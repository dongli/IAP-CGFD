#include "geomtk/Cartesian.h"

using namespace std;

int main(int argc, const char *argv[])
{
    ConfigManager configManager;
    Domain domain(1);
    Mesh mesh(domain);
    Field<double, 2> u, f;
    Field<double> fu, fstar, A, gamma, ustar;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;

    double dt, dx;
    const double eps = 1.0e-15;
    string outputPattern = "tspas.%3s.nc";

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }

    // Read configuration from file.
    configManager.parse(argv[1]);
    dt = configManager.getValue("tspas", "dt", 1.0);
    dx = configManager.getValue("tspas", "dx", 0.01);
    outputPattern = configManager.getValue("tspas", "output_pattern", outputPattern);

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/dx);

    // Set the time manager.
    Time startTime(Date(2000, 1, 1), Seconds(0));
    Time endTime(Date(2000, 1, 1), Seconds(200));
    timeManager.init(startTime, endTime, dt);

    // Set up velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 1, true);
    f.create("f", "kg m-1", "tracer density", mesh, CENTER, 1);
    fu.create("fu", "kg s-1", "tracer mass flux", mesh, X_FACE, 1);
    fstar.create("fstar", "kg m-1", "intermediate tracer density", mesh, CENTER, 1);
    ustar.create("ustar", "m s-1", "intermediate velocity component along x axis", mesh, X_FACE, 1);
    A.create("A", "kg2 m-2", "internal variable", mesh, CENTER, 1);
    gamma.create("gamma", "1", "internal variable", mesh, X_FACE, 1);

    // Set the initial conditions.
    newIdx = oldIdx+1;
    for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
        u(oldIdx, i) = 0.005;
        u(newIdx, i) = 0.005;
    }
    u.applyBndCond(oldIdx);
    u.applyBndCond(newIdx, true);
    for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
        const SpaceCoord &x = mesh.gridCoord(CENTER, i);
        if (x(0) >= 0.05 && x(0) <= 0.1) {
            f(oldIdx, i) = 1.0;
        } else {
            f(oldIdx, i) = 0.0;
        }
    }
    f.applyBndCond(oldIdx);

    // Set up IO manager.
    io.init(timeManager);
    outputFileIdx = io.addOutputFile(mesh, outputPattern, Seconds(dt));
    io.addField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&f});

    // Run the main loop.
    double C = dt/dx;
    while (!timeManager.isFinished()) {
        newIdx = oldIdx+1; halfIdx = oldIdx+0.5;
        // Lax-Wendroff pass.
        for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
            fu(i) = C*0.5*(      u(halfIdx, i)    *(f(oldIdx, i+1)+f(oldIdx, i))-
                           C*pow(u(halfIdx, i), 2)*(f(oldIdx, i+1)-f(oldIdx, i)));
        }
        fu.applyBndCond();
        // Calculate gamma.
        for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
            double alpha = fabs(u(halfIdx, i))*C;
            gamma(i) = alpha*(1-alpha);
        }
        gamma.applyBndCond();
        // Calculate intermediate fstar.
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            double beta = fmax(2/(2-gamma(i-1)), 2/(2-gamma(i)));
            fstar(i) = f(oldIdx, i)-beta*(fu(i)-fu(i-1));
        }
        fstar.applyBndCond();
        // Calculate A.
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            double fMax = fmax(fmax(f(oldIdx, i-1), f(oldIdx, i)), f(oldIdx, i+1));
            double fMin = fmin(fmin(f(oldIdx, i-1), f(oldIdx, i)), f(oldIdx, i+1));
            A(i) = (fstar(i)-fMax)*(fstar(i)-fMin);
        }
        A.applyBndCond();
        // Calculate ustar.
        for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
            double tmp1 = fabs(A(i  ))+A(i  );
            double tmp2 = fabs(A(i+1))+A(i+1);
            double cstar = 0.5 *(tmp1/(fabs(A(i))+eps)+tmp2/(fabs(A(i+1))+eps))-
                           0.25*(tmp1*tmp2/(fabs(A(i))*fabs(A(i+1))+eps));
            // NOTE: We let cstar be 0 or 1.
            cstar = cstar < 0.5 ? 0 : 1;
            ustar(i) = (cstar+(1-cstar)*C*fabs(u(halfIdx, i)))*u(halfIdx, i);
        }
        // Upwind pass.
        for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
            fu(i) = dt/dx*0.5*(u(halfIdx, i) *(f(oldIdx, i+1)+f(oldIdx, i))-
                               fabs(ustar(i))*(f(oldIdx, i+1)-f(oldIdx, i)));
        }
        fu.applyBndCond();
        // Get the final f.
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            f(newIdx, i) = f(oldIdx, i)-(fu(i)-fu(i-1));
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    return 0;
}
