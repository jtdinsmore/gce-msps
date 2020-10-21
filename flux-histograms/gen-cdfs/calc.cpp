#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <math.h>
#include <cmath>

typedef std::pair<double, double> CoordPair;
typedef std::pair<int, int> IndexPair;
typedef std::vector<std::vector<double>> DoubleVector;

//#define LOTS_OF_THREADS // Use on Erebus
/* When LOTS_OF_THREADS is defined, a target of 80-90 threads are made and run concurrently.
   Otherwise, four threads are made and run concurrently. */

#define SENSITIVITY_DIVISOR 20.0
#define SMOOTHING
#define SIGMA_TH 0.28
#define K_TH 0.45

#define pi 3.14159265358979323
#define ONE_PLUS_EPSILON 1.0000000001
#define max(a, b) (((a) > (b)) ? (a) : (b))

#define ALPHA 1.94
#define L_MIN 1.0e29
#define DIST_TO_CENTER 8.5 // kpc
#define CM_PER_KPC_SQUARED 9.523396e+42

#define CDF_BINS 2000
#define CDF_FLUX_LOW 1e-14
#define CDF_FLUX_HIGH 8e-11

#define ALPHA_HIST 1.94
#define NPTF_L_BREAK_HIST 2.5389429e+34
#define NPTF_N_1_HIST -0.66
#define NPTF_N_2_HIST 18.2
#define BARTELS_18_L_BREAK_HIST 1.7378008e+33
#define BARTELS_18_N_1_HIST 0.97
#define BARTELS_18_N_2_HIST 2.6



const double L_MIN_RANGE[2] = { 1.0e28, 1.0e34 };
const double L_MAX_RANGE[2] = {1.0e34, 1.0e38}; //{ 1.0e34, 1.0e36 };
const double ALPHA_RANGE[2] = { 1.1, 2.5 };
const double L0_RANGE[2] = {1.0e30, 2.0e36}; //{ 1.0e32, 2.0e34 };
const double SIGMA_RANGE[2] = { 0.001, 1 };

const double fermilabPaperPoint[2] = { 1e35, 1e29 };
const double bartels15Point[2] = { 1.5, 7.0e34 };
const double logNormalPaperPoint[2] = { 0.88e34, 0.62 };
const double logNormalPloegPoint[2] = { 1.3023e+32, 0.69550156 };
const double logNormalGautamPoint[2] = { 4.2970e+30, 0.93936155 };

// Values for the computation of the integral of the NFW profile.
#define GAMMA 1.2
#define NFW_SCALE_DIST 20 // kpc
#define INTEGRAL_WIDTH (20 * DIST_TO_CENTER)
#define NFW_INTEGRAL_STEPS 10000
#define PLOT_SIZE 50
#define FERMILAB_ALPHA 1.94

#define ROOT "/home/gridsan/jdinsmore/gce-gen-cdfs/"
#define SENSITIVITY_PATH (ROOT "sensitivity/sensitivity_10.txt")
#define PLOEG_PATH (ROOT "luminosity-models-step/ploeg/data/disk.csv")
#define DISPLAY_SIZE 20.0
#define myAbs(a) ((a) > 0 ? (a) : -(a))


enum class LUMINOSITY_FUNCTION {
    POWER_LAW,
    BARTELS_15,
    LOG_NORMAL,
    PLOEG,
    GAUTAM,
    NPTF,
    BARTELS_18,
};

class SensitivityMap {
public:
    SensitivityMap() {
        std::ifstream sensitivityFile;
        sensitivityFile.open(SENSITIVITY_PATH);
        if (sensitivityFile.is_open()) {
            std::string line;
            while (std::getline(sensitivityFile, line)) {
                std::vector<double> v;
                std::stringstream ss(line);
                for (double d; ss >> d;) {
                    v.push_back(d / SENSITIVITY_DIVISOR);
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
                fluxes.push_back(v);
            }
        }
        else {
            std::cout << "Sensitivity file not found." << std::endl;
        }

        upperLeft[0] = int((1 - DISPLAY_SIZE / 90.0) * fluxes.size() / 2.0);
        upperLeft[1] = int((1 - DISPLAY_SIZE / 180.0) * fluxes[0].size() / 2.0);
        lowerRight[0] = int((1 + DISPLAY_SIZE / 90.0) * fluxes.size() / 2.0);
        lowerRight[1] = int((1 + DISPLAY_SIZE / 180.0) * fluxes[0].size() / 2.0);
        skyShape = { lowerRight[0] - upperLeft[0], lowerRight[1] - upperLeft[1] };
    }

public:
    CoordPair indexToLatLon(IndexPair xy) const {
        double deltaXFrac = 1 - xy.first / (skyShape.first / 2.0);
        double deltaYFrac = xy.second / (skyShape.second / 2.0) - 1;
        return { DISPLAY_SIZE * deltaXFrac * pi / 180.0, DISPLAY_SIZE * deltaYFrac * pi / 180.0 };
    }
    IndexPair latLonToIndex(CoordPair latLon) const {
        double deltaXFrac = latLon.first * 180.0 / pi / DISPLAY_SIZE;
        double deltaYFrac = latLon.second * 180.0 / pi / DISPLAY_SIZE;
        return { int((1 - deltaXFrac) * (skyShape.first / 2.0)), int((1 + deltaYFrac) * (skyShape.second / 2.0)) };
    }
    double getFluxThreshold(CoordPair latLon) const {
        IndexPair xy = latLonToIndex(latLon);
        return fluxes[xy.first + upperLeft[0]][xy.second + upperLeft[1]];
    }

public:
    IndexPair skyShape;
    int upperLeft[2];
    int lowerRight[2];
    DoubleVector fluxes;
};

const SensitivityMap thresholds;

double lumFunc(LUMINOSITY_FUNCTION func, double lum) {
    switch(func) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
        return pow(lum, -ALPHA_HIST) * exp(-lum / fermilabPaperPoint[0]);
    case LUMINOSITY_FUNCTION::BARTELS_15:
        return pow(lum, -bartels15Point[0]) * exp(-lum / bartels15Point[1]);
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        return 1 / (lum) * exp(-pow(log10(lum) - log10(logNormalPaperPoint[0]), 2) / (2 * logNormalPaperPoint[1] * logNormalPaperPoint[1]));
    case LUMINOSITY_FUNCTION::PLOEG:
        return 1 / (lum) * exp(-pow(log10(lum) - log10(logNormalPloegPoint[0]), 2) / (2 * logNormalPloegPoint[1] * logNormalPloegPoint[1]));
    case LUMINOSITY_FUNCTION::GAUTAM:
        return 1 / (lum) * exp(-pow(log10(lum) - log10(logNormalGautamPoint[0]), 2) / (2 * logNormalGautamPoint[1] * logNormalGautamPoint[1]));
    case LUMINOSITY_FUNCTION::NPTF:
        if (lum < NPTF_L_BREAK_HIST) {
            return pow(lum / NPTF_L_BREAK_HIST, -NPTF_N_1_HIST);
        }
        else {
            return pow(lum / NPTF_L_BREAK_HIST, -NPTF_N_2_HIST);
        }
    case LUMINOSITY_FUNCTION::BARTELS_18:
        if (lum < BARTELS_18_L_BREAK_HIST) {
            return pow(lum / BARTELS_18_L_BREAK_HIST, -BARTELS_18_N_1_HIST);
        }
        else {
            return pow(lum / BARTELS_18_L_BREAK_HIST, -BARTELS_18_N_2_HIST);
        }
    default:
        throw "Invalid lum func";
    }
}

double sensitivity_func(double flux, double threshold) {
#ifndef SMOOTHING
    return (flux > threshold) ? 1 : 0;
#else
    return 0.5 * (1 + erf((log10(flux) - (log10(threshold) + K_TH)) / (sqrt(2) * SIGMA_TH)));
#endif
}

double pdfAtLatLon(LUMINOSITY_FUNCTION func, double flux, CoordPair latLon) {
    const double deltaRadialDistance = INTEGRAL_WIDTH / NFW_INTEGRAL_STEPS;
    double integral = 0;
    double distFromCenter, nfwSquaredValue, lum;
    const double cosLat = cos(latLon.first);
    const double cosLon = cos(latLon.second);
    double threshold = thresholds.getFluxThreshold(latLon);
    for (double radialDistance = deltaRadialDistance; radialDistance < INTEGRAL_WIDTH; radialDistance += deltaRadialDistance) {
        distFromCenter = sqrt(radialDistance * radialDistance + DIST_TO_CENTER * DIST_TO_CENTER
            - 2 * DIST_TO_CENTER * radialDistance * cosLon * cosLat);
        nfwSquaredValue = pow(distFromCenter / NFW_SCALE_DIST, -GAMMA) * pow(1 + distFromCenter / NFW_SCALE_DIST, -3 + GAMMA);
        nfwSquaredValue = nfwSquaredValue * nfwSquaredValue;

        lum = flux * 4 * pi * radialDistance * radialDistance * CM_PER_KPC_SQUARED;

        integral += sensitivity_func(flux, threshold) * nfwSquaredValue * deltaRadialDistance
            * radialDistance * radialDistance * radialDistance * radialDistance
            * cosLat * lumFunc(func, lum);
    }

    return integral;
}

double pdf(LUMINOSITY_FUNCTION func, double flux) {
    // returns an array of the value across the entire sky.
    // The values are scaled relative to each other in the image, but do not produce the correct total luminosity across the entire GCE
    double integral = 0;
    CoordPair latLon;
    for (int x = 0; x < thresholds.skyShape.first; x++) {
        //std::cout << x << '/' << skyMap->size() << std::endl;
        for (int y = 0; y < thresholds.skyShape.second; y++) {
            latLon = thresholds.indexToLatLon({ x, y });
            if (myAbs(latLon.first) > 2 * pi / 180) {// Cut out 2 degrees around the equator on each side.
                integral += pdfAtLatLon(func, flux, latLon);
            }
        }
    }

    return integral;
}

void writeCDF(LUMINOSITY_FUNCTION func, const double* fluxes) {
    double pdfs[CDF_BINS];
    double cdfs[CDF_BINS];
    for (int i = 0; i < CDF_BINS; i++) {
        std::cout << i << "/" << CDF_BINS << std::endl;
        pdfs[i] = pdf(func, fluxes[i]);
    }
    for (int i = 0; i < CDF_BINS; i++) {
        cdfs[i] = 0;
        for(int j = i; j < CDF_BINS; j++) {
            cdfs[i] += pdfs[j];
        }
    }

    std::ofstream f;
#ifdef SMOOTHING
    std::string tab = "smoothing";
#else
    std::string tab = "position";
#endif
    std::string name = "data-" + std::to_string(int(SENSITIVITY_DIVISOR)) +  "x/";
    switch(func) {
        case LUMINOSITY_FUNCTION::POWER_LAW:
            name += "power-law-" + tab + ".dat";
            std::cout << "Power law done" << std::endl;
            break;
        case LUMINOSITY_FUNCTION::BARTELS_15:
            name += "bartels15-" + tab + ".dat";
            std::cout << "Bartels 15 done" << std::endl;
            break;
        case LUMINOSITY_FUNCTION::LOG_NORMAL:
            name += "log-normal-" + tab + ".dat";
            std::cout << "Log normal done" << std::endl;
            break;
        case LUMINOSITY_FUNCTION::PLOEG:
            name += "ploeg-" + tab + ".dat";
            std::cout << "Ploeg done" << std::endl;
            break;
        case LUMINOSITY_FUNCTION::GAUTAM:
            name += "gautam-" + tab + ".dat";
            std::cout << "Gautam done" << std::endl;
            break;
        case LUMINOSITY_FUNCTION::NPTF:
            name += "nptf-" + tab + ".dat";
            std::cout << "NPTF done" << std::endl;
            break;
        case LUMINOSITY_FUNCTION::BARTELS_18:
            name += "bartels18-" + tab + ".dat";
            std::cout << "Bartels 18 done" << std::endl;
            break;
    }
    std::cout << name << std::endl;
    f.open(name);
    for (int i = 0; i < CDF_BINS; i++) {
        f << fluxes[i] << " " << pdfs[i] <<  " "  << cdfs[i] << std::endl;
    }
    f.close();
}

int main(int argc, char** argv) {
    std::cout << "Sensitivity " << SENSITIVITY_DIVISOR << std::endl;

    double fluxes[CDF_BINS];
    for (int i = 0; i < CDF_BINS; i++) {
        fluxes[i] = pow(10, log10(CDF_FLUX_LOW) + (double)(i) / CDF_BINS *
            (log10(CDF_FLUX_HIGH) - log10(CDF_FLUX_LOW)));
    }

    std::thread powerlaw = std::thread(writeCDF, LUMINOSITY_FUNCTION::POWER_LAW, &fluxes[0]);
    std::thread bartels15 = std::thread(writeCDF, LUMINOSITY_FUNCTION::BARTELS_15, &fluxes[0]);
    std::thread lognormal = std::thread(writeCDF, LUMINOSITY_FUNCTION::LOG_NORMAL, &fluxes[0]);
    std::thread ploeg = std::thread(writeCDF, LUMINOSITY_FUNCTION::PLOEG, &fluxes[0]);
    std::thread gautam = std::thread(writeCDF, LUMINOSITY_FUNCTION::GAUTAM, &fluxes[0]);
    std::thread nptf = std::thread(writeCDF, LUMINOSITY_FUNCTION::NPTF, &fluxes[0]);
    std::thread bartels18 = std::thread(writeCDF, LUMINOSITY_FUNCTION::BARTELS_18, &fluxes[0]);
    powerlaw.join();
    lognormal.join();
    ploeg.join();
    gautam.join();
    nptf.join();
    bartels15.join();
    bartels18.join();

    return 0;
}
