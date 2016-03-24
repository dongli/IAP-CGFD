#include "geomtk/Cartesian.h"

using namespace std;
using arma::field;

#define OMEGA 0.3

inline double mismatch(double fm1, double f, double fp1, const string &limiterType) {
    double df = (fp1-fm1)*0.5; // maybe 0.25?
    if (limiterType == "none") {
        return df;
    } else if (limiterType == "monotonic") {
        double dfMin = f-min(min(fm1, f), fp1);
        double dfMax = max(max(fm1, f), fp1)-f;
        return copysign(min(min(abs(df), dfMin), dfMax), df);
    } else if (limiterType == "positive_definite") {
        return copysign(min(abs(df), 2*f), df);
    } else {
        REPORT_ERROR("Unknown limiter type \"" << limiterType);
    }
}

inline void ppm(double fm2, double fm1, double f, double fp1, double fp2, const string &limiterType, double &df, double &f6) {
    double dfl = mismatch(fm2, fm1, f,   limiterType);
           df  = mismatch(fm1, f,   fp1, limiterType);
    double dfr = mismatch(f,   fp1, fp2, limiterType);
    double fl = 0.5*(fm1+f)+(dfl-df)/6.0;
    double fr = 0.5*(fp1+f)+(df-dfr)/6.0;
    fl = f-copysign(min(abs(df), abs(fl-f)), df);
    fr = f+copysign(min(abs(df), abs(fr-f)), df);
    f6 = 6*(f-0.5*(fl+fr));
    df = fr-fl;

    double zero = 0.5*df-0.25*f6;
    assert(zero == 0);
}

void callFluxOperator(const Mesh &mesh, double alpha, const string &fluxOperator,
                      const string &limiterType,
                      const field<double> &u, const field<double> &v,
                      const field<double> &fx, const field<double> &fy,
                      field<double> &fu, field<double> &fv) {
    field<double> dfx(size(fx)), fx6(size(fx)), dfy(size(fy)), fy6(size(fy));
    if (fluxOperator == "ppm") {
        // Calculate the subgrid distribution of tracer.
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                ppm(fx(i-2, j), fx(i-1, j), fx(i, j), fx(i+1, j), fx(i+2, j), limiterType, dfx(i, j), fx6(i, j));
                ppm(fy(i, j-2), fy(i, j-1), fy(i, j), fy(i, j+1), fy(i, j+2), limiterType, dfy(i, j), fy6(i, j));
            }
        }
    }
    // Along x axis.
    for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
            fu(i, j) = 0;
            double CN = u(i, j)*alpha;
            int K = static_cast<int>(CN), I = static_cast<int>(i+1-CN);
            double c = CN-K;
            // Calculate integer flux.
            // NOTE: The flux has direction.
            if (K >= 1) {
                for (int k = 1; k <= K; ++k) {
                    fu(i, j) += fx(i+1-k, j);
                }
            } else if (K <= -1) {
                for (int k = 1; k <= -K; ++k) {
                    fu(i, j) -= fx(i+k, j);
                }
            }
            // Calculate fractional flux.
            if (fluxOperator == "upwind") {
                fu(i, j) += c*fx(I, j);
            } else if (fluxOperator == "van_leer") {
                double df = mismatch(fx(I-1, j), fx(I, j), fx(I+1, j), limiterType);
                fu(i, j) += c*(fx(I, j)+(copysign(1, c)-c)*df*0.5);
            } else if (fluxOperator == "ppm") {
                double x1, x2;
                if (c > 0) {
                    x1 = 1-c; x2 = 1;
                } else {
                    x1 = 0; x2 = -c; I += 1;
                }
                double dx = x2-x1, dx2 = x2*x2-x1*x1, dx3 = x2*x2*x2-x1*x1*x1;
                fu(i, j) += copysign(fx(I, j)*dx+0.5*dfx(I, j)*dx2+fx6(I, j)*(dx/12-dx3/3), c);
            }
        }
    }
    // Along y axis.
    for (int j = mesh.js(HALF); j <= mesh.je(HALF); ++j) {
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            fv(i, j) = 0;
            double CN = v(i, j)*alpha;
            int K = static_cast<int>(CN), J = static_cast<int>(j+1-CN);
            double c = CN-K;
            // Calculate integer flux.
            // NOTE: The flux has direction.
            if (K >= 1) {
                for (int k = 1; k <= K; ++k) {
                    fv(i, j) += fy(i, j+1-k);
                }
            } else if (K <= -1) {
                for (int k = 1; k <= -K; ++k) {
                    fv(i, j) -= fy(i, j+k);
                }
            }
            // Calculate fractional flux.
            if (fluxOperator == "upwind") {
                fv(i, j) += c*fy(i, J);
            } else if (fluxOperator == "van_leer") {
                double df = mismatch(fy(i, J-1), fy(i, J), fy(i, J+1), limiterType);
                fv(i, j) += c*(fy(i, J)+(copysign(1, c)-c)*df*0.5);
            } else if (fluxOperator == "ppm") {
                double y1, y2;
                if (c > 0) {
                    y1 = 1-c; y2 = 1;
                } else {
                    y1 = 0; y2 = -c; J += 1;
                }
                double dy = y2-y1, dy2 = y2*y2-y1*y1, dy3 = y2*y2*y2-y1*y1*y1;
                fv(i, j) += copysign(fy(i, J)*dy+0.5*dfy(i, J)*dy2+fy6(i, J)*(dy/12-dy3/3), c);
            }
        }
    }
}

int main(int argc, char const* argv[])
{
    ConfigManager configManager;
    Domain domain(2);
    Mesh mesh(domain, 20);
    Field<double, 2> u, v, f;
    Field<double> fx, fy, cx, cy, au, av, fu, fv;
    TimeManager timeManager;
    IOManager io;
    Regrid regrid(mesh);
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;

    double dt, dx, dy;
    string outputPattern = "ffsl.%3s.nc";
    string advectiveOperator = "semi_lagrangian"; // semi_lagrangian, from_flux_operator
    string fluxOperator = "van_leer"; // upwind, van_leer, ppm
    string limiterType = "none"; // none, monotonic, positive_definite
    RegridMethod advectInterpMethod;

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }

    // Read configuration from file.
    configManager.parse(argv[1]);
    dt = configManager.getValue("ffsl", "dt", 0.05);
    dx = configManager.getValue("ffsl", "dx", 0.01);
    dy = configManager.getValue("ffsl", "dy", 0.01);
    outputPattern = configManager.getValue("ffsl", "output_pattern", outputPattern);
    advectiveOperator = configManager.getValue("ffsl", "advective_operator", advectiveOperator);
    fluxOperator = configManager.getValue("ffsl", "flux_operator", fluxOperator);
    limiterType = configManager.getValue("ffsl", "limiter_type", limiterType);
    advectInterpMethod = Regrid::methodFromString(configManager.getValue<string>("ffsl", "advect_interp_method"));

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m", 0, geomtk::BndType::PERIODIC, 1, geomtk::BndType::PERIODIC);
    domain.setAxis(1, "y", "y axis", "m", 0, geomtk::BndType::PERIODIC, 1, geomtk::BndType::PERIODIC);

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
    fx.create("fx", "kg m-2", "new density due to advective operator along x axis", mesh, CENTER, 2);
    fy.create("fy", "kg m-2", "new density due to advective operator along y axis", mesh, CENTER, 2);
    cx.create("cx", "1", "divergence component along x axis times dt", mesh, CENTER, 2);
    cy.create("cy", "1", "divergence component along x axis times dt", mesh, CENTER, 2);
    au.create("au", "kg m-2", "advective operator along x axis", mesh, X_FACE, 2);
    av.create("av", "kg m-2", "advective operator along y axis", mesh, Y_FACE, 2);
    fu.create("fu", "kg m-2", "flux operator along x axis", mesh, X_FACE, 2);
    fv.create("fv", "kg m-2", "flux operator along y axis", mesh, Y_FACE, 2);

    // Set initial conditions.
    // NOTE: Velocity does not satisfy periodic boundary condition!
    SpaceCoord x0(2);
    x0.set(0.25, 0.5);
    newIdx = oldIdx+1;
    for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            const SpaceCoord &x = mesh.gridCoord(CENTER, i, j);
            double d = domain.calcDistance(x, x0);
            u(oldIdx, i, j) = -OMEGA*(x(1)-0.5);
            v(oldIdx, i, j) =  OMEGA*(x(0)-0.5);
            u(oldIdx, i, j) = 0.134;
            v(oldIdx, i, j) = 0;
            u(newIdx, i, j) = u(oldIdx, i, j);
            v(newIdx, i, j) = v(oldIdx, i, j);
            if (fabs(x(0)-x0(0)) >= 0.02 && d < 0.1) {
                f(oldIdx, i, j) = 1;
            } else if (fabs(x(0)-x0(0)) <= 0.02 && d < 0.1 && x(1) >= 0.55) {
                f(oldIdx, i, j) = 1;
            } else {
                f(oldIdx, i, j) = 0;
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
        if (advectiveOperator == "semi_lagrangian") {
            // Use a low-order semi-Lagrangian method to calculate the intermediate values.
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    const SpaceCoord &x0 = mesh.gridCoord(CENTER, i, j);
                    SpaceCoord xd(2);
                    MeshIndex meshIdx(2);
                    double fx0, fy0;
                    // Along x axis.
                    double ua = 0.5*(u(halfIdx, i, j)+u(halfIdx, i-1, j));
                    xd(0) = x0(0)-dt*ua; xd(1) = x0(1); // first-order trajectory
                    domain.validateCoord(xd);
                    meshIdx.locate(mesh, xd);
                    regrid.run(advectInterpMethod, oldIdx, f, xd, fx0, &meshIdx);
                    fx(i, j) = 0.5*(f(oldIdx, i, j)+fx0);
                    // Along y axis.
                    double va = 0.5*(v(halfIdx, i, j)+v(halfIdx, i, j-1));
                    xd(1) = x0(1)-dt*va; xd(0) = x0(0);
                    domain.validateCoord(xd);
                    meshIdx.locate(mesh, xd);
                    regrid.run(advectInterpMethod, oldIdx, f, xd, fy0, &meshIdx);
                    fy(i, j) = 0.5*(f(oldIdx, i, j)+fy0);
                }
            }
        } else if (advectiveOperator == "from_flux_operator") {
            // Use the advective form of the outer flux operator.
            callFluxOperator(mesh, dt/dx, fluxOperator, limiterType, u(halfIdx), v(halfIdx), f(oldIdx), f(oldIdx), au(), av());
            // Substract the divergence parts from the flux.
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    cx(i, j) = dt/dx*(u(halfIdx, i, j)-u(halfIdx, i-1, j));
                    cy(i, j) = dt/dy*(v(halfIdx, i, j)-v(halfIdx, i, j-1));
                }
            }
            cx.applyBndCond(); cy.applyBndCond();
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
                    au(i, j) += 0.5*(cx(i, j)*f(oldIdx, i, j)+cx(i+1, j)*f(oldIdx, i+1, j));
                }
            }
            for (int j = mesh.js(HALF); j <= mesh.je(HALF); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    av(i, j) += 0.5*(cy(i, j)*f(oldIdx, i, j)+cy(i, j+1)*f(oldIdx, i, j+1));
                }
            }
            au.applyBndCond(); av.applyBndCond();
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    fx(i, j) = f(oldIdx, i, j)-0.5*(au(i, j)-au(i-1, j));
                    fy(i, j) = f(oldIdx, i, j)-0.5*(av(i, j)-av(i, j-1));
                }
            }
        }
        fx.applyBndCond(); fy.applyBndCond();
        // Do outer flux operator.
        callFluxOperator(mesh, dt/dx, fluxOperator, limiterType, u(halfIdx), v(halfIdx), fy(), fx(), fu(), fv());
        fu.applyBndCond(); fv.applyBndCond();
        // Get the final results.
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                f(newIdx, i, j) = f(oldIdx, i, j)-(fu(i, j)-fu(i-1, j))-(fv(i, j)-fv(i, j-1));
            }
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    return 0;
}
