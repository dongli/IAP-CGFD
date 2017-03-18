#include "geomtk/Cartesian.h"

using namespace std;
using arma::field;

inline double mismatch(double fm1, double f, double fp1, const string &limiterType) {
    double df = (fp1-fm1)*0.5;
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

inline void ppm(double fm2, double fm1, double f, double fp1, double fp2, const string &limiterType, double &fl, double &df, double &f6) {
    double dfl = mismatch(fm2, fm1, f,   limiterType);
           df  = mismatch(fm1, f,   fp1, limiterType);
    double dfr = mismatch(f,   fp1, fp2, limiterType);
           fl = 0.5*(fm1+f)+(dfl-df)/6.0;
    double fr = 0.5*(fp1+f)+(df-dfr)/6.0;
    fl = f-copysign(min(abs(df), abs(fl-f)), df);
    fr = f+copysign(min(abs(df), abs(fr-f)), df);
    f6 = 6*f-3*(fl+fr);
    df = fr-fl;
}

void callFluxOperator(const Mesh &mesh, double alpha, const string &fluxOperator,
                      const string &limiterType, const field<double> &u,
                      const field<double> &f, 
                      Field<double> &fl, Field<double> &df,
                      Field<double> &f6, Field<double> &fu) {
    if (fluxOperator == "ppm") {
        // Calculate the subgrid distribution of tracer.
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            ppm(f(i-2), f(i-1), f(i), f(i+1), f(i+2), limiterType, fl(i), df(i), f6(i));
        }
        fl.applyBndCond();
        df.applyBndCond();
        f6.applyBndCond();
    }
    for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
        fu(i) = 0;
        double CN = u(i)*alpha;
        int K = static_cast<int>(CN);
        double c = CN-K;
        // Calculate integer flux.
        // NOTE: The flux has direction.
        if (K >= 1) {
            for (int k = 1; k <= K; ++k) {
                fu(i) += f(i+1-k);
            }
        } else if (K <= -1) {
            for (int k = 1; k <= -K; ++k) {
                fu(i) -= f(i+k);
            }
        }
        int I = CN > 0 ? i-K : i+1-K;
        // Calculate fractional flux.
        if (fluxOperator == "upwind") {
            fu(i) += c*f(I);
        } else if (fluxOperator == "van_leer") {
            double df = mismatch(f(I-1), f(I), f(I+1), limiterType);
            fu(i) += c*(f(I)+(copysign(1, c)-c)*df*0.5);
        } else if (fluxOperator == "ppm") {
            double x1, x2;
            if (c > 0) {
                x1 = 1-c; x2 = 1;
            } else {
                x1 = 0; x2 = -c;
            }
            double dx = x2-x1, dx2 = x2*x2-x1*x1, dx3 = x2*x2*x2-x1*x1*x1;
            fu(i) += copysign(fl(I)*dx+0.5*df(I)*dx2+f6(I)*(0.5*dx2-dx3/3.0), c); // This formula is from original Colella and Woodward (1984).
            // fu(i) += copysign(f(I)*dx+0.5*df(I)*dx2+f6(I)*(dx/12-dx3/3), c); // This formula in Carpenter et al. (1990) is wrong!
        }
    }
    fu.applyBndCond();
}

int main(int argc, char const* argv[])
{
    ConfigManager configManager;
    Domain domain(1);
    Mesh mesh(domain, 2);
    Field<double, 2> u, f;
    Field<double> fl, df, f6, fu;
    TimeManager timeManager;
    IOManager io;
    Regrid regrid(mesh);
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;

    double dt, dx;
    string outputPattern = "ffsl.%3s.nc";
    string fluxOperator = "van_leer"; // upwind, van_leer, ppm
    string limiterType = "none"; // none, monotonic, positive_definite

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }

    // Read configuration from file.
    configManager.parse(argv[1]);
    dt = configManager.getValue("ffsl_1d", "dt", 1.00);
    dx = configManager.getValue("ffsl_1d", "dx", 0.01);
    outputPattern = configManager.getValue("ffsl_1d", "output_pattern", outputPattern);
    fluxOperator = configManager.getValue("ffsl_1d", "flux_operator", fluxOperator);
    limiterType = configManager.getValue("ffsl_1d", "limiter_type", limiterType);

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m", 0, geomtk::BndType::PERIODIC, 1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/dx);

    // Set the time manager.
    Time startTime(Date(2000, 1, 1), Seconds(0));
    Time endTime(Date(2000, 1, 1), Seconds(200));
    timeManager.init(startTime, endTime, dt);

    // Create velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 1, true);
    f.create("f", "kg m-2", "tracer density", mesh, CENTER, 1);
    fl.create("fl", "", "tracer density on left cell face", mesh, CENTER, 1);
    df.create("df", "", "tracer density mismatch", mesh, CENTER, 1);
    f6.create("f6", "", "curvature", mesh, CENTER, 1);
    fu.create("fu", "kg m-2", "flux operator along x axis", mesh, X_FACE, 1);

    // Set initial conditions.
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

    // Set IO manager.
    io.init(timeManager);
    outputFileIdx = io.addOutputFile(mesh, outputPattern, Seconds(dt));
    io.addField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&f});

    // Run the main loop.
    while (!timeManager.isFinished()) {
        halfIdx = oldIdx+0.5; newIdx = oldIdx+1;
        callFluxOperator(mesh, dt/dx, fluxOperator, limiterType, u(halfIdx), f(oldIdx), fl, df, f6, fu);
        // Get the final results.
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            f(newIdx, i) = f(oldIdx, i)-(fu(i)-fu(i-1));
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    return 0;
}
