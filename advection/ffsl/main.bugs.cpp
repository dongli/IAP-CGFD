#include "geomtk/Cartesian.h"

using namespace std;

#define OMEGA 0.3

int main(int argc, char const* argv[])
{
    ConfigManager configManager;
    Domain domain(2);
    Mesh mesh(domain, 20);
    Field<double, 2> u, v, f;
    Field<double> fx, fy, fu, fv;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;
    
    Regrid regrid(mesh);

    double dt, dx, dy;
    string outputPattern = "ffsl.%3s.nc";
    string advectiveOperator = "semi_lagrangian";
    string fluxOperator = "van_leer"; // upwind, van_leer, ppm
    RegridMethod advectInterpMethod;

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }

    // Read configuration from file.
    configManager.parse(argv[1]);
    dt = configManager.getValue("ffsl", "dt", 0.5);
    dx = configManager.getValue("ffsl", "dx", 0.01);
    dy = configManager.getValue("ffsl", "dy", 0.01);
    outputPattern = configManager.getValue("ffsl", "output_pattern", outputPattern);
    advectiveOperator = configManager.getValue("ffsl", "advective_operator", advectiveOperator);
    fluxOperator = configManager.getValue("ffsl", "flux_operator", fluxOperator);
    advectInterpMethod = Regrid::methodFromString(configManager.getValue<string>("ffsl", "advect_interp_method"));

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);
    domain.setAxis(1, "y", "y axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/dx, domain.axisSpan(1)/dy);

    // Set the time manager.
    Time startTime(Date(2000, 1, 1), Seconds(0));
    Time endTime(Date(2000, 1, 1), Seconds(6*PI2/OMEGA));
    timeManager.init(startTime, endTime, dt);

    // Create velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 2, true);
    v.create("v", "m s-1", "velocity component along y axis", mesh, Y_FACE, 2, true);
    f.create("f", "kg m-2", "tracer density", mesh, CENTER, 2);
    fx.create("fx", "kg m-2", "new tracer density due to x advective inner operator", mesh, CENTER, 2);
    fy.create("fy", "kg m-2", "new tracer density due to y advective inner operator", mesh, CENTER, 2);
    fu.create("fu", "kg s-1", "tracer mass flux along x axis", mesh, X_FACE, 2);
    fv.create("fv", "kg s-1", "tracer mass flux along y axis", mesh, Y_FACE, 2);

    // Set initial conditions.
    // NOTE: Velocity does not satisfy periodic boundary condition!
    SpaceCoord x0(2);
    x0.set(0.25, 0.5);
    newIdx = oldIdx+1;
    for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            const SpaceCoord &x = mesh.gridCoord(CENTER, i, j);
            double d = domain.calcDistance(x, x0);
            u(oldIdx)(i, j) = -OMEGA*(x(1)-0.5);
            v(oldIdx)(i, j) =  OMEGA*(x(0)-0.5);
            u(newIdx)(i, j) = u(oldIdx)(i, j);
            v(newIdx)(i, j) = v(oldIdx)(i, j);
            if (fabs(x(0)-x0(0)) >= 0.02 && d < 0.1) {
                f(oldIdx)(i, j) = 1;
            } else if (fabs(x(0)-x0(0)) <= 0.02 && d < 0.1 && x(1) >= 0.55) {
                f(oldIdx)(i, j) = 1;
            } else {
                f(oldIdx)(i, j) = 0;
            }
        }
    }
    u.applyBndCond(oldIdx); u.applyBndCond(newIdx, true);
    v.applyBndCond(oldIdx); v.applyBndCond(newIdx, true);
    f.applyBndCond(oldIdx);

    // Set IO manager.
    io.init(timeManager);
    outputFileIdx = io.addOutputFile(mesh, outputPattern, Seconds(dt));
    io.addField(outputFileIdx, "double", FULL_DIMENSION, {&u, &v, &f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&u, &v, &f});

    // Run the main loop.
    while (!timeManager.isFinished()) {
        halfIdx = oldIdx+0.5; newIdx = oldIdx+1;
        // Do inner advective operator.
        // Use a low-order 1D advective-form semi-Lagrangian scheme.
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                const SpaceCoord &x0 = mesh.gridCoord(CENTER, i, j);
                SpaceCoord xd(2);
                MeshIndex meshIdx(2);
                double fx0, fy0;
                if (advectiveOperator == "semi_lagrangian") {
                    // Along x axis.
                    double ua = 0.5*(u(halfIdx)(i, j)+u(halfIdx)(i-1, j));
                    xd(0) = x0(0)-dt*ua; xd(1) = x0(1); // first-order trajectory
                    domain.validateCoord(xd);
                    meshIdx.locate(mesh, xd);
                    regrid.run(advectInterpMethod, oldIdx, f, xd, fx0, &meshIdx);
                    fx()(i, j) = 0.5*(f(oldIdx)(i, j)+fx0);
                    // Along y axis.
                    double va = 0.5*(v(halfIdx)(i, j)+v(halfIdx)(i, j-1));
                    xd(1) = x0(1)-dt*va; xd(0) = x0(0);
                    domain.validateCoord(xd);
                    meshIdx.locate(mesh, xd);
                    regrid.run(advectInterpMethod, oldIdx, f, xd, fy0, &meshIdx);
                    fy()(i, j) = 0.5*(f(oldIdx)(i, j)+fy0);
                }
            }
        }
        fx.applyBndCond(); fy.applyBndCond();
        // Do outer flux operator.
        // Along x axis.
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
                fu()(i, j) = 0;
                double CN = u(halfIdx)(i, j)*dt/dx;
                int K = static_cast<int>(CN), I = i+1-K;
                double c = CN-K;
                // Calculate integer flux.
                if (K >= 1) {
                    for (int k = 1; k <= K; ++k) {
                        fu()(i, j) += fy()(i+1-k, j);
                    }
                } else if (K <= -1) {
                    for (int k = 1; k <= -K; ++k) {
                        fu()(i, j) += fy()(i+k, j);
                    }
                }
                // Calculate fractional flux.
                if (fluxOperator == "van_leer") {
                    fu()(i, j) += c*(fy()(I, j)+0.25*(fy()(I+1, j)-fy()(I-1, j))*(copysign(1, c)-c));
                }
            }
        }
        // Along y axis.
        for (int j = mesh.js(HALF); j <= mesh.je(HALF); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                fv()(i, j) = 0;
                double CN = v(halfIdx)(i, j)*dt/dy;
                int K = static_cast<int>(CN), J = j+1-K;
                double c = CN-K;
                // Calculate integer flux.
                if (K >= 1) {
                    for (int k = 1; k <= K; ++k) {
                        fv()(i, j) += fx()(i, j+1-k);
                    }
                } else if (K <= -1) {
                    for (int k = 1; k <= -K; ++k) {
                        fv()(i, j) += fx()(i, j+k);
                    }
                }
                // Calculate fractional flux.
                if (fluxOperator == "van_leer") {
                    fv()(i, j) += c*(fx()(i, J)+0.25*(fx()(i, J+1)-fx()(i, J-1))*(copysign(1, c)-c));
                }
            }
        }
        // Get the final results.
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                f(newIdx)(i, j) = f(oldIdx)(i, j)-(fu()(i, j)-fu()(i-1, j))-(fv()(i, j)-fv()(i, j-1));
            }
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    return 0;
}
