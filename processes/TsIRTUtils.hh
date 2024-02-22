#ifndef TsIRTUtils_hh
#define TsIRTUtils_hh

#include "TsIRTConfiguration.hh"
#include "globals.hh"
#include <vector>
#include <map>

class TsIRTUtils {
public:
    TsIRTUtils();
   ~TsIRTUtils();

    // IRT Sampling Helper Functions
    G4double NormQuantile(G4double x);
    G4double erfcx_y100(G4double x);
    G4double erfcx(G4double x);
    G4double erfc(G4double x);
    G4double erfcInv(G4double x);
    G4double erfcWxy(G4double c, G4double x, G4double y);
    G4double SampleTypeII(G4double alpha, G4double sigma, G4double r0, G4double D);
    G4double SamplePDC(G4double a, G4double b);
    G4double Lambda(G4double x, G4double beta, G4double alphatilde);
    
    // Bin Index Finding Functions
    G4int FindBin(G4double x, std::vector<G4double> vx);
    G4int FindBin(G4int n, G4double xmin, G4double xmax, G4double);

    // Output Time Recording
    std::vector<G4double> CreateTimeSteps(G4double t0, G4double tf, G4int tbins, G4bool isLogSpace);

    // Temperature Related Functions
    G4double WaterDielectricConstant(G4double TKelvin);
    G4double OnsagerRadius(G4double TKelvin);
    G4double DebyeFactor(G4double TKelvin, G4double ChargeA, G4double ChargeB, G4double Radius);
    G4double WaterDensity(G4double TCelsius);
    G4double NoyesRelationship(G4double Kobs, G4double Kact, G4double Kdiff);
    G4double ArrheniusFunction(G4double A, G4double T, G4double E);
    G4double SmoluchowskiFunction(G4double Beta, G4double Radius, G4double Diffusion);
    G4double DebyeFunction(G4double Beta, G4double Radius, G4double Diffusion, G4double TCelsius, G4double ChargeA, G4double ChargeB);

    // Temperature Dependent Diffusion Rates
    G4double ElectronDiffusionRate(G4double TCelsius);
    G4double H3ODiffusionRate(G4double TKelvin);
    G4double OHmDiffusionRate(G4double TKelvin);
    G4double WaterDiffusionRate(G4double TKelvin, G4double Diff);

    // Polynomial Root Solver
    G4int roots(double*,int, double*, double*);
    void deflate(double*,int, double*,double*, double *);
    void find_quad(double*, int, double*, double*, double*, int*);
    void diff_poly(double *, int, double*);
    void recurse(double *, int, double *, int, double*, double *, int *);
    void get_quads(double*, int, double*, double*);
    std::vector<double> GetRoots(int, std::vector<double>);
    
};

#endif

