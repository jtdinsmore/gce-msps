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
#include <boost/math/special_functions/gamma.hpp>

#define INVALID -3789213.4812047 // A random number designed to encode an invalid input.

typedef std::pair<double, double> CoordPair;
typedef std::pair<int, int> IndexPair;
typedef std::vector<std::vector<double>> DoubleVector;

//#define LOTS_OF_THREADS // Use on Erebus
/* When LOTS_OF_THREADS is defined, a target of 80-90 threads are made and run concurrently.
   Otherwise, four threads are made and run concurrently. */

#define SENSITIVITY_DIVISOR 1.0
#define DATA_RELEASE 2

#define CUT_UNRESOLVED
#define CUT_RESOLVED
#define CUT_RADIUS 2.0

#define pi 3.14159265358979323
#define ONE_PLUS_EPSILON 1.0000000001
#define max(a, b) (((a) > (b)) ? (a) : (b))

#define ALPHA 1.94
#define L_MIN 1.0e29
#define DIST_TO_CENTER 8.5 // kpc
#define CM_PER_KPC_SQUARED 9.523396e+42
#ifndef CUT_UNRESOLVED
    #define FLUX_EXCESS 1.794113925439598e-09  // All flux units are in ergs per second per square centimeter
    #define NPTF_L_BREAK 2.5389429e+34
#else
    #define FLUX_EXCESS 1.27508891111732e-09
    #define NPTF_L_BREAK 2.4127028e+34

#endif

#define BARTELS_18_L_BREAK 1.7378008e+33


#define REAL_HISTOGRAM

#ifdef REAL_HISTOGRAM
    #define NUM_STEPS_WITHIN_BIN 100
    #define FLUX_HIST_BINS 7
    #define FLUX_HIST_LOW 1.9590460265411667e-12
    #define FLUX_HIST_HIGH  1.637482454933416e-11
#else
    #define NUM_STEPS_WITHIN_BIN 1
    #define FLUX_HIST_BINS 100
    #define FLUX_HIST_LOW 1.0e-13
    #define FLUX_HIST_HIGH 8.0e-11
#endif
#define ALPHA_HIST 1.94
#define BARTELS_15_ALPHA_HIST 1.5
#define NPTF_N_1_HIST -0.66
#define NPTF_N_2_HIST 18.2


#define BARTELS_18_N_1_HIST 0.97
#define BARTELS_18_N_2_HIST 2.60
#define BARTELS_18_L_MIN 1.0e30
#define BARTELS_18_L_MAX 1.0e37

#define BARTELS_15_L_MIN 1.0e29





std::array<double, FLUX_HIST_BINS> fluxValues, fluxWidths, fluxLow, fluxHigh;

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
#define NFW_INTEGRAL_STEPS 1000 // 10000 FOR CALCULATIONS, 1000 FOR GRAPHS.
#define PLOT_SIZE 50
#define FERMILAB_ALPHA 1.94

#define ROOT "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/"
#if DATA_RELEASE == 1
    #define SENSITIVITY_PATH (ROOT "sensitivity/sensitivity_8_mask.txt")
#elif DATA_RELEASE == 2 
    #define SENSITIVITY_PATH (ROOT "sensitivity/sensitivity_10.txt")
#endif
#define PLOEG_PATH (ROOT "luminosity-models-step/ploeg/data/disk.csv")
#define DISPLAY_SIZE 20.0
#define myAbs(a) ((a) > 0 ? (a) : -(a))

double nptfLBreak;
double nptfLMin = 1e29;

enum class VALUE {
    TOTAL_NUM,
    TOTAL_FLUX,
    SEEN_NUM,
    SEEN_FLUX,
};

enum class LUMINOSITY_FUNCTION {
    POWER_LAW,
    BARTELS_15,
    LOG_NORMAL,
    PLOEG,
    GAUTAM,
    NPTF,
    BARTELS_18,
    POWER_LAW_ALPHA,
    ERROR,
};

enum class METHOD {
    COUNT,
    HIST,
    ERROR,
};

LUMINOSITY_FUNCTION luminosityFunction = LUMINOSITY_FUNCTION::POWER_LAW;

std::string sciNot(double value) {
    int power = log10(value);
    if (log10(value) < 0 && value != pow(10, power)) {
        power--;
    }
    if (myAbs(power) > 6) {
        return std::to_string(value / pow(10.0, power)) + "e" + std::to_string(power);
    }
    return std::to_string(value);
}

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

class PloegData {
public:
    PloegData() {
        std::ifstream ploegFile;
        ploegFile.open(PLOEG_PATH);
        if (ploegFile.is_open()) {
            std::string line;
            while (std::getline(ploegFile, line)) {
                std::vector<double> v;
                std::stringstream ss(line);
                for (double d; ss >> d;) {
                    v.push_back(d); // L, dN/dlogL
                    if (ss.peek() == ',') {
                        ss.ignore();
                    }
                }
                xs.push_back(pow(10, v[0]));
                logxs.push_back(v[0]);
                ys.push_back(log10(exp(1)) / pow(10, v[0]) * v[1]);
            }
        }
        else {
            std::cout << "Ploeg data not found." << std::endl;
        }

        normalization = integrate(pow(10, logxs[0]) * ONE_PLUS_EPSILON, false);
    }
    double operator()(double l) const {// l is the log10 of the luminosity
        if (l <= logxs[0]) { return ys[0]; }
        if (l >= logxs[logxs.size() - 1]) { return ys[ys.size() - 1]; }

        // Binary search for l:
        int leftIndex = 0;
        int rightIndex = logxs.size() - 1;
        for (;;) {
            int midIndex = (rightIndex + leftIndex) / 2;
            if (logxs[midIndex] > l) {
                rightIndex = midIndex;
            }
            else if (logxs[midIndex] < l) {
                leftIndex = midIndex;
            }
            else {
                return ys[midIndex];
            }
            if (rightIndex - leftIndex <= 1) {
                assert(leftIndex != rightIndex);
                return ys[leftIndex] + (ys[rightIndex] - ys[leftIndex]) * (l - logxs[leftIndex]) / (logxs[rightIndex] - logxs[leftIndex]);
            }
        }
    }

public:
    double integrate(double start, bool normalize=true) const {
        int i;
        for (i = 0; i < xs.size() && xs[i] < start; i++);
        assert(i < xs.size());
        assert(i > 0);
        double frac = (start - xs[i - 1]) / (xs[i] - xs[i - 1]);
        double startY = ys[i - 1] + (ys[i] - ys[i - 1]) * frac;
        double integral = 0.5 * (startY + ys[i]) * (xs[i] - start);
        for (; i < xs.size() - 1; i++) {
            integral += 0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i]);
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }

        /*double integral = 0;
        for (double logx = log10(start); logx < logxs[logxs.size() - 1]; logx += INTEGRATION_LOG_STEP) {
            // Trapezoidal integration
            double xLow = pow(10, logx);
            double xHigh = pow(10, logx + INTEGRATION_LOG_STEP);
            double yLow = operator()(logx);
            double yHigh = operator()(logx + INTEGRATION_LOG_STEP);
            integral += (yLow + yHigh) / 2 * (xHigh - xLow);
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }*/
    }

    double lintegrate(double start, bool normalize = true) const {int i;
        for (i = 0; i < xs.size() && xs[i] < start; i++);
        assert(i < xs.size());
        assert(i > 0);
        double frac = (start - xs[i - 1]) / (xs[i] - xs[i - 1]);
        double startY = ys[i - 1] + (ys[i] - ys[i - 1]) * frac;
        double integral = 1 / 6.0 * (-start * start * (ys[i] + 2 * startY) + start * xs[i] * (startY - ys[i]) + xs[i] * xs[i] * (2 * ys[i] + startY));
        for (; i < xs.size() - 1; i++) {
            integral += 1 / 6.0 * (ys[i] * (-2 * xs[i] * xs[i] + xs[i] * xs[i + 1] + xs[i + 1] * xs[i + 1])
                + ys[i + 1] * (-xs[i] * xs[i] - xs[i] * xs[i + 1] + 2 * xs[i + 1] * xs[i + 1]));
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }
        /*double integral = 0;
        for (double logx = log10(start); logx < logxs[logxs.size() - 1]; logx += INTEGRATION_LOG_STEP) {
            // Trapezoidal integration
            double xLow = pow(10, logx);
            double xHigh = pow(10, logx + INTEGRATION_LOG_STEP);
            double yLow = operator()(logx);
            double yHigh = operator()(logx + INTEGRATION_LOG_STEP);
            integral += (xLow * yLow + xHigh * yHigh) / 2 * (xHigh - xLow);
        }
        if (normalize) {
            return integral / normalization;
        }
        else {
            return integral;
        }*/
    }

public:
    std::vector<double> logxs, ys, xs;
    double normalization;
};

const SensitivityMap thresholds;
//const PloegData ploegData;

// ================= These functions characterize the luminosity function ============

// The luminosity function must be normalized.

double improvedGamma(double s, double x) {
    assert(x >= 0);
    while (s < 0) {
        return (improvedGamma(s+1, x) - pow(x, s) * exp(-x)) / s;
    }
    return boost::math::tgamma(s, x);
}

double integrate(double start, double arg1, double arg2) {// arg1, arg2; l0, sigma; nBelow, nAbove
    const double norm = (-arg2 * nptfLBreak + pow(BARTELS_18_L_MAX, 1 - arg2) * pow(nptfLBreak, arg2) -
        BARTELS_18_L_MIN * pow(nptfLBreak/BARTELS_18_L_MIN, arg1) + arg2 * BARTELS_18_L_MIN * pow(nptfLBreak/BARTELS_18_L_MIN, arg1) +
        arg1 * (nptfLBreak - pow(BARTELS_18_L_MAX, 1 - arg2) * pow(nptfLBreak, arg2)))/((-1 +
        arg1) * (-1 + arg2));

    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
        return improvedGamma(1 - ALPHA, start / arg2) / improvedGamma(1 - ALPHA, arg1 / arg2);
    case LUMINOSITY_FUNCTION::BARTELS_15:
        return improvedGamma(1 - BARTELS_15_ALPHA_HIST, start / arg2) / improvedGamma(1 - BARTELS_15_ALPHA_HIST, arg1 / arg2);
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
    case LUMINOSITY_FUNCTION::PLOEG:
    case LUMINOSITY_FUNCTION::GAUTAM:
        return 0.5 * (1 - erf(log10(start / arg1)/(sqrt(2) * arg2)));
    //case LUMINOSITY_FUNCTION::PLOEG:
    //    return ploegData.integrate(start);
    case LUMINOSITY_FUNCTION::NPTF:
        if (start < nptfLBreak) {
            return 1 + (-1 + arg2) * pow(nptfLBreak / start, arg1-1) / (arg1 - arg2);
        }
        else {
            return (-1 + arg1) * pow(nptfLBreak / start, arg2-1) / (arg1 - arg2);
        }
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        return improvedGamma(1 - arg1, start / arg2) / improvedGamma(1 - arg1, L_MIN / arg2);
    case LUMINOSITY_FUNCTION::BARTELS_18:
        if (start == INVALID) { start = BARTELS_18_L_MIN; }
        if (start > BARTELS_18_L_MAX) {return 0; }
        if (start < nptfLBreak) {
            return 1 / norm * (-arg2 * nptfLBreak + pow(BARTELS_18_L_MAX, 1 - arg2) * pow(nptfLBreak, arg2) -
                start * pow(nptfLBreak/start, arg1) + arg2 * start * pow(nptfLBreak/start, arg1) +
                arg1 * (nptfLBreak - pow(BARTELS_18_L_MAX, 1 - arg2) * pow(nptfLBreak, arg2)))/((-1 +
                arg1) * (-1 + arg2));
        }
        else {
            return 1/norm * (pow(BARTELS_18_L_MAX * start, -arg2) * (pow(BARTELS_18_L_MAX, arg2) * start -
               BARTELS_18_L_MAX * pow(start, arg2)) * pow(nptfLBreak, arg2))/(-1 + arg2);
        }
    default:
        std::cout << "Running integrate with an invalid luminosity function." << std::endl;
        return 0;
    }
}

double lintegrate(double start, double arg1, double arg2) {// lMin, lMax; l0, sigma
    const double norm = (-arg2 * nptfLBreak + pow(BARTELS_18_L_MAX, 1 - arg2) * pow(nptfLBreak, arg2) -
        BARTELS_18_L_MIN * pow(nptfLBreak/BARTELS_18_L_MIN, arg1) + arg2 * BARTELS_18_L_MIN * pow(nptfLBreak/BARTELS_18_L_MIN, arg1) +
        arg1 * (nptfLBreak - pow(BARTELS_18_L_MAX, 1 - arg2) * pow(nptfLBreak, arg2)))/((-1 +
        arg1) * (-1 + arg2));

    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
        if (start == INVALID) { start = arg1; }
        return arg2 * improvedGamma(2 - ALPHA, start / arg2) / improvedGamma(1 - ALPHA, arg1 / arg2);
    case LUMINOSITY_FUNCTION::BARTELS_15:
        if (start == INVALID) { start = arg1; }
        return arg2 * improvedGamma(2 - BARTELS_15_ALPHA_HIST, start / arg2) / improvedGamma(1 - BARTELS_15_ALPHA_HIST, arg1 / arg2);
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
    case LUMINOSITY_FUNCTION::PLOEG:
    case LUMINOSITY_FUNCTION::GAUTAM:
        if (start == INVALID) { start = 0; }
        return 0.5 * arg1 * exp(arg2 * arg2 * log(10)*log(10) / 2) * (1 - erf((log10(start) - log10(arg1) - arg2 * arg2 * log(10)) / (sqrt(2) * arg2)));
    //case LUMINOSITY_FUNCTION::PLOEG:
    //    if (start == INVALID) { start = pow(10, ploegData.logxs[0]) * ONE_PLUS_EPSILON; }
    //    return ploegData.lintegrate(start);
    case LUMINOSITY_FUNCTION::NPTF:
        if (start == INVALID) { start = 0; }
        if (start < nptfLBreak) {
            return  (arg1 - arg2) / (arg1 - arg2 - (1 - arg2) * pow(nptfLMin / nptfLBreak, 1 - arg1))
                * nptfLBreak * (1 - arg1) * (1-arg2) * (1 / ((arg1-2) * (arg2-2)) + pow(nptfLBreak / start, arg1 - 2) / ((arg1 - 2) * (arg1 - arg2)));
        }
        else {
            return  (arg1 - arg2) / (arg1 - arg2 - (1 - arg2) * pow(nptfLMin / nptfLBreak, 1 - arg1))
                * nptfLBreak * (1 - arg1) * (1 - arg2) * (pow(nptfLBreak / start, arg2 - 2) / ((arg2 - 2) * (arg1 - arg2)));
        }
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        if (start == INVALID) { start = L_MIN; }
        return arg2 * improvedGamma(2 - arg1, start / arg2) / improvedGamma(1 - arg1, L_MIN / arg2);
    case LUMINOSITY_FUNCTION::BARTELS_18:
        if (start == INVALID) { start = BARTELS_18_L_MIN; }
        if (start > BARTELS_18_L_MAX) {return 0; }
        if (start < nptfLBreak) {
            return 1 / norm * (-arg2 * pow(nptfLBreak, 2) - 2 * pow(start, 2 - arg1) * pow(nptfLBreak, arg1) +
                arg2 * pow(start, 2 - arg1) * pow(nptfLBreak, arg1) +
                2 * pow(BARTELS_18_L_MAX, 2 - arg2) * pow(nptfLBreak, arg2) +
                arg1 * (pow(nptfLBreak, 2) - pow(BARTELS_18_L_MAX, 2 - arg2) * pow(nptfLBreak, arg2)))/((-2 +
                arg1) * (-2 + arg2));
        }
        else {
            return 1/norm * ((pow(BARTELS_18_L_MAX, 2 - arg2) - pow(start, 2 - arg2)) * pow(nptfLBreak, arg2))/(2 - arg2);
        }
    default:
        std::cout << "Running lintegrate with an invalid luminosity function." << std::endl;
        return 0;
    }
}


// ================= Declare helper functions =======================

double numSeenFunc(double threshold, double arg1, double arg2) {// Returns(unscaled) number of visible pulsars
    //assert(threshold > arg1);
    return integrate(threshold, arg1, arg2);
}

double fluxSeenFunc(double threshold, double arg1, double arg2) {// Returns(unscaled) amount of luminosity visible
    //assert(threshold > arg1);
    return lintegrate(threshold, arg1, arg2);
}

double totalFluxFunc(double threshold, double arg1, double arg2) {// Returns(unscaled) total luminosity
    return lintegrate(INVALID, arg1, arg2);
}


double getValueAtLatLon(CoordPair latLon, double fluxThreshold, VALUE value, double arg1, double arg2) {
    // returns the (unscaled) value for a specific lon and lat
    // Integrate value along the line of sight:
    const double deltaRadialDistance = INTEGRAL_WIDTH / NFW_INTEGRAL_STEPS;
    double integral = 0;
    double distFromCenter, nfwSquaredValue, lThreshold;
    const double cosLat = cos(latLon.first);
    const double cosLon = cos(latLon.second);
    for (double radialDistance = deltaRadialDistance; radialDistance < INTEGRAL_WIDTH; radialDistance += deltaRadialDistance) {
        distFromCenter = sqrt(radialDistance * radialDistance + DIST_TO_CENTER * DIST_TO_CENTER
            - 2 * DIST_TO_CENTER * radialDistance * cosLon * cosLat);

        #ifdef CUT_UNRESOLVED
        if ((value == VALUE::TOTAL_FLUX || value == VALUE::TOTAL_NUM) && distFromCenter > CUT_RADIUS) {continue;}
        #endif

        #ifdef CUT_RESOLVED
        if ((value == VALUE::SEEN_FLUX || value == VALUE::SEEN_NUM) && distFromCenter > CUT_RADIUS) {continue;}
        #endif

        nfwSquaredValue = pow(distFromCenter / NFW_SCALE_DIST, -GAMMA) * pow(1 + distFromCenter / NFW_SCALE_DIST, -3 + GAMMA);
        nfwSquaredValue = nfwSquaredValue * nfwSquaredValue;

        lThreshold = fluxThreshold * 4 * pi * radialDistance * radialDistance * CM_PER_KPC_SQUARED;

        switch (value) {
        case VALUE::TOTAL_NUM:
            integral += nfwSquaredValue * deltaRadialDistance * radialDistance * radialDistance * cosLat * CM_PER_KPC_SQUARED;
            break;
        case VALUE::SEEN_NUM:
            integral += nfwSquaredValue * deltaRadialDistance * radialDistance * radialDistance * cosLat * numSeenFunc(lThreshold, arg1, arg2) * CM_PER_KPC_SQUARED;
            break;
        case VALUE::TOTAL_FLUX:
            integral += nfwSquaredValue * deltaRadialDistance * cosLat / (4 * pi) * totalFluxFunc(lThreshold, arg1, arg2);
            break;
        case VALUE::SEEN_FLUX:
            integral += nfwSquaredValue * deltaRadialDistance * cosLat / (4 * pi) * fluxSeenFunc(lThreshold, arg1, arg2);
            break;
        }
    }

    return integral;
}

double get_lum_value(double lum) {
    double prob;
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        prob = pow(lum, -ALPHA_HIST) * exp(-lum / fermilabPaperPoint[0]) / (improvedGamma(1 - ALPHA_HIST, fermilabPaperPoint[1] / fermilabPaperPoint[0]) * pow(fermilabPaperPoint[0], 1 - ALPHA_HIST));
        if (lum < fermilabPaperPoint[1]) {
            prob = 0;
        }
        break;
    case LUMINOSITY_FUNCTION::BARTELS_15:
        prob = pow(lum, -BARTELS_15_ALPHA_HIST) * exp(-lum / bartels15Point[1]) / (improvedGamma(1 - BARTELS_15_ALPHA_HIST, BARTELS_15_L_MIN / bartels15Point[1]) * pow(bartels15Point[1], 1 - BARTELS_15_ALPHA_HIST));
        if (lum < BARTELS_15_L_MIN) {
            prob = 0;
        }
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        prob = log10(exp(1)) / (logNormalPaperPoint[1] * sqrt(2 * pi) * lum) * exp(-pow(log10(lum) - log10(logNormalPaperPoint[0]), 2) / (2 * logNormalPaperPoint[1] * logNormalPaperPoint[1]));
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        prob = log10(exp(1)) / (logNormalPloegPoint[1] * sqrt(2 * pi) * lum) * exp(-pow(log10(lum) - log10(logNormalPloegPoint[0]), 2) / (2 * logNormalPloegPoint[1] * logNormalPloegPoint[1]));
        break;
    case LUMINOSITY_FUNCTION::GAUTAM:
        prob = log10(exp(1)) / (logNormalGautamPoint[1] * sqrt(2 * pi) * lum) * exp(-pow(log10(lum) - log10(logNormalGautamPoint[0]), 2) / (2 * logNormalGautamPoint[1] * logNormalGautamPoint[1]));
        break;

    case LUMINOSITY_FUNCTION::NPTF:
        if (lum < NPTF_L_BREAK) {
            prob = pow(lum / NPTF_L_BREAK, -NPTF_N_1_HIST);
        }
        else {
            prob = pow(lum / NPTF_L_BREAK, -NPTF_N_2_HIST);
        }
        prob *= (1 - NPTF_N_1_HIST) * (1 - NPTF_N_2_HIST) / (NPTF_L_BREAK * (NPTF_N_1_HIST - NPTF_N_2_HIST));
        break;

    case LUMINOSITY_FUNCTION::BARTELS_18:
        if (lum < BARTELS_18_L_MIN || lum > BARTELS_18_L_MAX) {
            prob = 0;
            break;
        }

        if (lum < BARTELS_18_L_BREAK) {
            prob = pow(lum / BARTELS_18_L_BREAK, -BARTELS_18_N_1_HIST);
        }
        else {
            prob = pow(lum / BARTELS_18_L_BREAK, -BARTELS_18_N_2_HIST);
        }
        prob /= (-BARTELS_18_N_2_HIST * BARTELS_18_L_BREAK + pow(BARTELS_18_L_MAX, 1 - BARTELS_18_N_2_HIST) * pow(BARTELS_18_L_BREAK, BARTELS_18_N_2_HIST) -
            BARTELS_18_L_MIN * pow(BARTELS_18_L_BREAK/BARTELS_18_L_MIN, BARTELS_18_N_1_HIST) + BARTELS_18_N_2_HIST * BARTELS_18_L_MIN * pow(BARTELS_18_L_BREAK/BARTELS_18_L_MIN, BARTELS_18_N_1_HIST) +
            BARTELS_18_N_1_HIST * (BARTELS_18_L_BREAK - pow(BARTELS_18_L_MAX, 1 - BARTELS_18_N_2_HIST) * pow(BARTELS_18_L_BREAK, BARTELS_18_N_2_HIST)))/((-1 +
            BARTELS_18_N_1_HIST) * (-1 + BARTELS_18_N_2_HIST));

        break;
    }
    return prob;
}

std::array<double, FLUX_HIST_BINS> getHistAtLatLon(CoordPair latLon, double fluxThreshold) {
    // returns the (unscaled) value for a specific lon and lat

    // Integrate value along the line of sight:
    const double deltaRadialDistance = INTEGRAL_WIDTH / NFW_INTEGRAL_STEPS;
    std::array<double, FLUX_HIST_BINS> integral;
    integral.fill(0);
    double distFromCenter, nfwSquaredValue, lThreshold;
    const double cosLat = cos(latLon.first);
    const double cosLon = cos(latLon.second);

    for (double radialDistance = deltaRadialDistance; radialDistance < INTEGRAL_WIDTH; radialDistance += deltaRadialDistance) {
        distFromCenter = sqrt(radialDistance * radialDistance + DIST_TO_CENTER * DIST_TO_CENTER
            - 2 * DIST_TO_CENTER * radialDistance * cosLon * cosLat);
        nfwSquaredValue = pow(distFromCenter / NFW_SCALE_DIST, -GAMMA) * pow(1 + distFromCenter / NFW_SCALE_DIST, -3 + GAMMA);
        nfwSquaredValue = nfwSquaredValue * nfwSquaredValue;
        double fluxToLum = 4 * pi * radialDistance * radialDistance * CM_PER_KPC_SQUARED;

        #ifdef CUT_RESOLVED
        if (distFromCenter > CUT_RADIUS) {continue;}
        #endif

        lThreshold = fluxThreshold * fluxToLum;

        for (int i = 0; i < FLUX_HIST_BINS; i++) {
            double small_width = fluxWidths[i] * fluxToLum / NUM_STEPS_WITHIN_BIN;
            double prob = 0;
            for (int j = 0; j < NUM_STEPS_WITHIN_BIN; j++) {
                double lum = (fluxLow[i] + fluxWidths[i] * j / NUM_STEPS_WITHIN_BIN) * fluxToLum;
                if (lum > lThreshold) {
                    prob += get_lum_value(lum) * small_width;
                }
            }
            integral[i] += nfwSquaredValue * radialDistance * radialDistance * deltaRadialDistance * prob * cosLat; // Treating nfw here as a number density
        }
    }
    return integral;
}

void generateValueSkyMap(VALUE value, DoubleVector* skyMap, double arg1, double arg2) {
    // returns an array of the value across the entire sky.
    // The values are scaled relative to each other in the image, but do not produce the correct total luminosity across the entire GCE
    *skyMap = DoubleVector(thresholds.skyShape.first, std::vector<double>(thresholds.skyShape.second, 0));
    for (int x = 0; x < skyMap->size(); x++) {
        for (int y = 0; y < (*skyMap)[x].size(); y++) {
            CoordPair latLon = thresholds.indexToLatLon({ x, y });
            if (myAbs(latLon.first) > 2 * pi / 180) {// Cut out 2 degrees around the equator on each side.
                double val = getValueAtLatLon(latLon, thresholds.getFluxThreshold(latLon), value, arg1, arg2);
                (*skyMap)[x][y] = val;
            }
            else {
                (*skyMap)[x][y] = 0;
            }
        }
    }
}

double getValueAtConfig(VALUE value, double arg1, double arg2) {
    // Return the(unscaled) value at a certain arg1 or arg2
    DoubleVector skyMap;
    generateValueSkyMap(value, &skyMap, arg1, arg2);
    double sum = 0;
    for (int x = 0; x < skyMap.size(); x++) {
        for (int y = 0; y < skyMap[x].size(); y++) {
            sum += skyMap[x][y];// The values in the sky map have already been multiplied by the position distro.
        }
    }
    return sum;
}

void generatePowerLawPlotMap(VALUE value, DoubleVector* plotMap) {
    // Generate unscaled data
    const CoordPair powerSteps = { pow(L_MIN_RANGE[1] / L_MIN_RANGE[0], 1.0 / PLOT_SIZE),
                pow(L_MAX_RANGE[1] / L_MAX_RANGE[0], 1.0 / PLOT_SIZE) };

#ifndef LOTS_OF_THREADS
    for (int i = 0; i < PLOT_SIZE; i++) {
        double lMin = L_MIN_RANGE[0] * pow(powerSteps.first, i);
        std::cout << "LMax: " << i << std::endl;    
        for (int j = 0; j < PLOT_SIZE; j++) {
            double lMax = L_MAX_RANGE[0] * pow(powerSteps.second, j);
            (*plotMap)[i][j] = getValueAtConfig(value, lMin, lMax);
        }
    }

#else
    for (int i = 0; i < PLOT_SIZE; i++) {
        double lMin = L_MIN_RANGE[0] * pow(powerSteps.first, i);
        std::vector<std::future<double>*> futures(PLOT_SIZE, nullptr);
        for (int j = 0; j < PLOT_SIZE; j++) {
            double lMax = L_MAX_RANGE[0] * pow(powerSteps.second, j);
            //futures[j] = new std::future<double>(std::move(std::async(std::launch::async, &getValueAtConfig, value, lMin, lMax)));
        }
        std::cout << PLOT_SIZE << " threads made for value " << (int)value << "." << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            (*plotMap)[i][j] = 1; //futures[j]->get();
            //delete futures[j];
        }
    }
    std::cout << "Value " << (int)value << " completed." << std::endl;
#endif
}

void generatePowerLawAlphaPlotMap(VALUE value, DoubleVector* plotMap) {
    //Return the(unscaled) plot map
    CoordPair powerSteps = { (ALPHA_RANGE[1] - ALPHA_RANGE[0]) / PLOT_SIZE,
                pow(L_MAX_RANGE[1] / L_MAX_RANGE[0], 1.0 / PLOT_SIZE) };

    *plotMap = DoubleVector(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));

#ifndef LOTS_OF_THREADS
    for (int i = 0; i < PLOT_SIZE; i++) {
        double alpha = i * powerSteps.first + ALPHA_RANGE[0];
        //std::cout << alpha << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            double lMax = L_MAX_RANGE[0] * pow(powerSteps.second, j);

            (*plotMap)[i][j] = getValueAtConfig(value, alpha, lMax);
        }
    }

#else
    std::cout << "Unimplemented" << std::endk;
    throw "Unimplemented";
#endif
}

void generateLogNormalPlotMap(VALUE value, DoubleVector* plotMap) {
    //Return the(unscaled) plot map
    const CoordPair powerSteps = { pow(L0_RANGE[1] / L0_RANGE[0], 1.0 / PLOT_SIZE),
                    (SIGMA_RANGE[1] - SIGMA_RANGE[0]) / PLOT_SIZE };
    *plotMap = DoubleVector(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));

#ifndef LOTS_OF_THREADS
    for (int i = 0; i < PLOT_SIZE; i++) {
        double l0 = L0_RANGE[0] * pow(powerSteps.first, i);
        //std::cout << i << " / " << plotShape.first << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            double sigma = SIGMA_RANGE[0] + j * powerSteps.second;
            (*plotMap)[i][j] = getValueAtConfig(value, l0, sigma);
        }
    }

#else
    for (int i = 0; i < PLOT_SIZE; i++) {
        double l0 = L0_RANGE[0] * pow(powerSteps.first, i);
        std::vector<std::future<double>*> futures(PLOT_SIZE, nullptr);
        for (int j = 0; j < PLOT_SIZE; j++) {
            double sigma = SIGMA_RANGE[0] + powerSteps.second * j;
            //futures[j] = new std::future<double>(std::move(std::async(std::launch::async, &getValueAtConfig, value, l0, sigma)));
        }
        std::cout << PLOT_SIZE << " threads made for value " << (int)value << "at l0 " << l0 <<"." << std::endl;
        for (int j = 0; j < PLOT_SIZE; j++) {
            (*plotMap)[i][j] = 1;//futures[j]->get();
            //delete futures[j];
        }
    }
    std::cout << "Value " << (int)value << " completed." << std::endl;
#endif
}

// ========================= Generate data =========================
int powerLaw() {
    std::cout << "Using a power law luminosity function" << std::endl << std::endl;

    std::string writeText, bartels15Text;

    double paperScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, fermilabPaperPoint[1], fermilabPaperPoint[0]);
    writeText = "Paper values: (lMin = " + sciNot(fermilabPaperPoint[1]) + " ergs/s, lMax = " + sciNot(fermilabPaperPoint[0]) + " ergs/s)\n" +
        "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, fermilabPaperPoint[1], fermilabPaperPoint[0]) * paperScale) +
        "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, fermilabPaperPoint[1], fermilabPaperPoint[0]) * paperScale) +
        "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, fermilabPaperPoint[1], fermilabPaperPoint[0]) * paperScale / FLUX_EXCESS) + "\n";
    std::cout << paperScale << std::endl;
    std::cout << writeText << std::endl;
    std::ofstream recordFile;


    luminosityFunction = LUMINOSITY_FUNCTION::POWER_LAW_ALPHA;

    double bartels15Scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, bartels15Point[0], bartels15Point[1]);
    bartels15Text = "Bartels 15 values: (alpha = " + sciNot(bartels15Point[0]) + ", lMax = " + sciNot(bartels15Point[1]) + " ergs/s)\n" +
        "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, bartels15Point[0], bartels15Point[1]) * bartels15Scale) +
        "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, bartels15Point[0], bartels15Point[1]) * bartels15Scale) +
        "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, bartels15Point[0], bartels15Point[1]) * bartels15Scale / FLUX_EXCESS) + "\n";
    std::cout << bartels15Scale << std::endl;
    std::cout << bartels15Text << std::endl;

    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/record-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    recordFile << writeText;
    recordFile << bartels15Text;
    recordFile.close();
    

    // Generate unscaled data
    const CoordPair powerSteps = { pow(L_MIN_RANGE[1] / L_MIN_RANGE[0], 1.0 / PLOT_SIZE),
                pow(L_MAX_RANGE[1] / L_MAX_RANGE[0], 1.0 / PLOT_SIZE) };
    DoubleVector totalNum(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));
    DoubleVector numSeen(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));
    DoubleVector fluxSeen(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));
    DoubleVector totalFlux(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));

#ifndef LOTS_OF_THREADS
    std::thread totalNumThread(generatePowerLawPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawPlotMap, VALUE::SEEN_NUM, &numSeen);
    std::thread lumSeenThread(generatePowerLawPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    totalNumThread.join();
    numSeenThread.join();
    lumSeenThread.join();
    totalLumThread.join();
#else
    // Do two at a time, because 50 * 2 = 100, and Erebus uses 80 threads at a time./
    // Of course, the number of threads is weakly less than 100.
    std::thread totalNumThread(generatePowerLawPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawPlotMap, VALUE::SEEN_NUM, &numSeen);
    totalNumThread.join();
    numSeenThread.join();
    std::thread lumSeenThread(generatePowerLawPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    lumSeenThread.join();
    totalLumThread.join();
#endif

    // Write data to an output file
    std::ofstream totalNumFile, numSeenFile, lumSeenFile;
    totalNumFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/total-num-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    numSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/num-seen-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    lumSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/lum-seen-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    for (int x = 0; x < totalNum.size(); x++) {
        for (int y = 0; y < totalNum[x].size(); y++) {
            double scale = FLUX_EXCESS / totalFlux[x][y];
            totalNumFile << totalNum[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            numSeenFile << numSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            lumSeenFile << fluxSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
        }
        totalNumFile << std::endl;
        numSeenFile << std::endl;
        lumSeenFile << std::endl;
    }

    std::cout << "Done." << std::endl;
    system("python \"" ROOT "luminosity-models-position/graph-power-law.py\"");
    return 1;
}

int powerLawAlpha() {
    std::cout << "Using a power law alpha luminosity function" << std::endl << std::endl;
    double paperScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, FERMILAB_ALPHA, fermilabPaperPoint[0]);
    std::string writeText = "Paper values: (alpha = " + std::to_string((int)FERMILAB_ALPHA) + ", lMax = " + sciNot(fermilabPaperPoint[0]) + " ergs/s)\n" +
        "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, FERMILAB_ALPHA, fermilabPaperPoint[0]) * paperScale) +
        "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, FERMILAB_ALPHA, fermilabPaperPoint[0]) * paperScale) +
        "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, FERMILAB_ALPHA, fermilabPaperPoint[0]) * paperScale / FLUX_EXCESS) + "\n";

    std::cout << writeText << std::endl;
    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/record-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    recordFile << writeText;

    // Generate unscaled data// Generate unscaled data
    const CoordPair powerSteps = { pow(L_MIN_RANGE[1] / L_MIN_RANGE[0], 1.0 / PLOT_SIZE),
                pow(L_MAX_RANGE[1] / L_MAX_RANGE[0], 1.0 / PLOT_SIZE) };
    DoubleVector totalNum(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));
    DoubleVector numSeen(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));
    DoubleVector fluxSeen(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));
    DoubleVector totalFlux(PLOT_SIZE, std::vector<double>(PLOT_SIZE, 0));

#ifndef LOTS_OF_THREADS
    std::thread totalNumThread(generatePowerLawAlphaPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawAlphaPlotMap, VALUE::SEEN_NUM, &numSeen);
    std::thread lumSeenThread(generatePowerLawAlphaPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawAlphaPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    totalNumThread.join();
    numSeenThread.join();
    lumSeenThread.join();
    totalLumThread.join();
#else
    // Do two at a time, because 50 * 2 = 100, and Erebus uses 80 threads at a time./
    // Of course, the number of threads is weakly less than 100.
    std::thread totalNumThread(generatePowerLawPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generatePowerLawPlotMap, VALUE::SEEN_NUM, &numSeen);
    totalNumThread.join();
    numSeenThread.join();
    std::thread lumSeenThread(generatePowerLawPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generatePowerLawPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    lumSeenThread.join();
    totalLumThread.join();
#endif

    // Write data to an output file
    std::ofstream totalNumFile, numSeenFile, lumSeenFile;
    totalNumFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/total-num-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    numSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/num-seen-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    lumSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law-alpha/lum-seen-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    for (int x = 0; x < totalNum.size(); x++) {
        for (int y = 0; y < totalNum[x].size(); y++) {
            double scale = FLUX_EXCESS / totalFlux[x][y];
            totalNumFile << totalNum[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            numSeenFile << numSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            lumSeenFile << fluxSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
        }
        totalNumFile << std::endl;
        numSeenFile << std::endl;
        lumSeenFile << std::endl;
    }

    std::cout << "Done." << std::endl;
    system("python \"" ROOT "luminosity-models-position/graph-power-law-alpha.py\"");
    return 1;
}

int logNormal() {
    std::cout << "Using a log normal luminosity function" << std::endl << std::endl;
    double paperScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalPaperPoint[0], logNormalPaperPoint[1]);;
    std::string paperText = "Paper values: (l0 = " + sciNot(logNormalPaperPoint[0]) + " ergs/s, sigma = " + sciNot(logNormalPaperPoint[1]) + ")\n"
        + "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, logNormalPaperPoint[0], logNormalPaperPoint[1]) * paperScale)
        + "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, logNormalPaperPoint[0], logNormalPaperPoint[1]) * paperScale)
        + "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, logNormalPaperPoint[0], logNormalPaperPoint[1]) * paperScale / FLUX_EXCESS) + "\n";

    double ploegScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalPloegPoint[0], logNormalPloegPoint[1]);
    std::string ploegText = "Ploeg values: (l0 = " + sciNot(logNormalPloegPoint[0]) + " ergs/s, sigma = " + sciNot(logNormalPloegPoint[1]) + ")\n"
        + "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, logNormalPloegPoint[0], logNormalPloegPoint[1]) * ploegScale)
        + "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, logNormalPloegPoint[0], logNormalPloegPoint[1]) * ploegScale)
        + "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, logNormalPloegPoint[0], logNormalPloegPoint[1]) * ploegScale / FLUX_EXCESS) + "\n";

    double gautamScale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalGautamPoint[0], logNormalGautamPoint[1]);
    std::string gautamText = "Gautam values: (l0 = " + sciNot(logNormalGautamPoint[0]) + " ergs/s, sigma = " + sciNot(logNormalGautamPoint[1]) + ")\n"
        + "\tTotal number of pulsars: " + sciNot(getValueAtConfig(VALUE::TOTAL_NUM, logNormalGautamPoint[0], logNormalGautamPoint[1]) * gautamScale)
        + "\n\tNumber of visible pulsars: " + sciNot(getValueAtConfig(VALUE::SEEN_NUM, logNormalGautamPoint[0], logNormalGautamPoint[1]) * gautamScale)
        + "\n\tFraction of seen luminosity: " + sciNot(getValueAtConfig(VALUE::SEEN_FLUX, logNormalGautamPoint[0], logNormalGautamPoint[1]) * gautamScale / FLUX_EXCESS) + "\n";

    std::cout << paperText << std::endl;
    std::cout << ploegText << std::endl;
    std::cout << gautamText << std::endl;
    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/record-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    recordFile << paperText << std::endl;
    recordFile << ploegText << std::endl;
    recordFile << gautamText << std::endl;

    // Generate unscaled data
    DoubleVector totalNum, numSeen, fluxSeen, totalFlux;

#ifndef LOTS_OF_THREADS
    std::thread totalNumThread(generateLogNormalPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generateLogNormalPlotMap, VALUE::SEEN_NUM, &numSeen);
    std::thread lumSeenThread(generateLogNormalPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generateLogNormalPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    totalNumThread.join();
    numSeenThread.join();
    lumSeenThread.join();
    totalLumThread.join();
#else
    // Do two at a time, because 50 * 2 = 100, and Erebus uses 80 threads at a time./
    // Of course, the number of threads is weakly less than 98.
    std::thread totalNumThread(generateLogNormalPlotMap, VALUE::TOTAL_NUM, &totalNum);
    std::thread numSeenThread(generateLogNormalPlotMap, VALUE::SEEN_NUM, &numSeen);
    totalNumThread.join();
    numSeenThread.join();
    std::thread lumSeenThread(generateLogNormalPlotMap, VALUE::SEEN_FLUX, &fluxSeen);
    std::thread totalLumThread(generateLogNormalPlotMap, VALUE::TOTAL_FLUX, &totalFlux);
    lumSeenThread.join();
    totalLumThread.join();
#endif

    // Write data to an output file
    std::ofstream totalNumFile, numSeenFile, lumSeenFile;
    totalNumFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/total-num-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    numSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/num-seen-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    lumSeenFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/lum-seen-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    for (int x = 0; x < totalNum.size(); x++) {
        for (int y = 0; y < totalNum[x].size(); y++) {
            double scale = FLUX_EXCESS / totalFlux[x][y];
            totalNumFile << totalNum[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            numSeenFile << numSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
            lumSeenFile << fluxSeen[x][y] * scale << (y == totalNum[x].size() - 1 ? "" : ", ");
        }
        totalNumFile << std::endl;
        numSeenFile << std::endl;
        lumSeenFile << std::endl;
    }

    std::cout << "Done." << std::endl;
    system("python \"" ROOT "luminosity-models-position/graph-log-normal.py\"");
    return 1;
}

int ploeg() {
    std::cout << "Using Ploeg et al.'s luminosity function" << std::endl << std::endl;

    std::future<double> totalFlux = std::async(std::launch::async, &getValueAtConfig, VALUE::TOTAL_FLUX, 0, 0);
    std::future<double> totalNum = std::async(std::launch::async, &getValueAtConfig, VALUE::TOTAL_NUM, 0, 0);
    std::future<double> seenFlux = std::async(std::launch::async, &getValueAtConfig, VALUE::SEEN_FLUX, 0, 0);
    std::future<double> seenNum = std::async(std::launch::async, &getValueAtConfig, VALUE::SEEN_NUM, 0, 0);

    double ploegScale = FLUX_EXCESS / totalFlux.get();
    double totalNumScaled = totalNum.get() * ploegScale;
    double numVisibleScaled = seenNum.get() * ploegScale;
    double fracFluxScaled = seenFlux.get() * ploegScale / FLUX_EXCESS;

    std::string ploegText = "Ploeg luminosity function \n"
        "\tTotal number of pulsars: " + sciNot(totalNumScaled)
        + "\n\tNumber of visible pulsars: " + sciNot(numVisibleScaled)
        + "\n\tFraction of seen luminosity: " + sciNot(fracFluxScaled) + "\n";

    std::cout << ploegText << std::endl;
    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/ploeg/record-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    recordFile << ploegText << std::endl;

    return 1;
}

int nptf() {
    std::cout << "Using the NPTF luminosity function" << std::endl << std::endl;
    std::future<double> totalFlux, totalNum, seenFlux, seenNum;
    double nptfScale, totalNumScaled, numVisibleScaled, fracFluxScaled;

    nptfLBreak = NPTF_L_BREAK;

    totalFlux = std::async(std::launch::async, getValueAtConfig, VALUE::TOTAL_FLUX, -0.66, 18.2);
    totalNum = std::async(std::launch::async, getValueAtConfig, VALUE::TOTAL_NUM, -0.66, 18.2);
    seenFlux = std::async(std::launch::async, getValueAtConfig, VALUE::SEEN_FLUX, -0.66, 18.2);
    seenNum = std::async(std::launch::async, getValueAtConfig, VALUE::SEEN_NUM, -0.66, 18.2);

    nptfScale = FLUX_EXCESS / totalFlux.get();
    totalNumScaled = totalNum.get() * nptfScale;
    numVisibleScaled = seenNum.get() * nptfScale;
    fracFluxScaled = seenFlux.get() * nptfScale / FLUX_EXCESS;

    std::string nptfText = "NPTF NFW PS luminosity function (nBelow = -0.66, nAbove = 18.2, lumBreak = " + sciNot(nptfLBreak) + " ergs/s \n"
        "\tTotal number of pulsars: " + sciNot(totalNumScaled)
        + "\n\tNumber of visible pulsars: " + sciNot(numVisibleScaled)
        + "\n\tFraction of seen luminosity: " + sciNot(fracFluxScaled) + "\n\n";
    std::cout << nptfText << std::endl;


    luminosityFunction = LUMINOSITY_FUNCTION::BARTELS_18;
    nptfLBreak = BARTELS_18_L_BREAK;

    std::cout << integrate(BARTELS_18_L_MAX/2, 0.97, 2.60) << std::endl;
    std::cout << lintegrate(BARTELS_18_L_MAX/2, 0.97, 2.60) << std::endl;

    totalFlux = std::async(std::launch::async, getValueAtConfig, VALUE::TOTAL_FLUX, 0.97, 2.60);
    totalNum = std::async(std::launch::async, getValueAtConfig, VALUE::TOTAL_NUM, 0.97, 2.60);
    seenFlux = std::async(std::launch::async, getValueAtConfig, VALUE::SEEN_FLUX, 0.97, 2.60);
    seenNum = std::async(std::launch::async, getValueAtConfig, VALUE::SEEN_NUM, 0.97, 2.60);

    nptfScale = FLUX_EXCESS / totalFlux.get();
    totalNumScaled = totalNum.get() * nptfScale;// Positive
    numVisibleScaled = seenNum.get() * nptfScale;// Negative
    fracFluxScaled = seenFlux.get() * nptfScale / FLUX_EXCESS;// Negative

    std::cout << nptfScale << std::endl; // Positive

    std::string bartelsText = "Bartels luminosity function (nBelow = 0.97, nAbove = 2.60, lumBreak = " + sciNot(nptfLBreak) + " ergs/s \n"
        "\tTotal number of pulsars: " + sciNot(totalNumScaled)
        + "\n\tNumber of visible pulsars: " + sciNot(numVisibleScaled)
        + "\n\tFraction of seen luminosity: " + sciNot(fracFluxScaled) + "\n\n";
    std::cout << bartelsText << std::endl;


    std::ofstream recordFile;
    recordFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/nptf/record-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
    //recordFile << nptfText << std::endl;
    recordFile << nptfText << bartelsText << std::endl;

    return 1;
}


int hist() {
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        std::cout << "Using a power law luminosity function" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        std::cout << "Using a log normal luminosity function" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        std::cout << "Using a log normal luminosity function (Ploeg)" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::GAUTAM:
        std::cout << "Using a log normal luminosity function (Gautam)" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::NPTF:
        std::cout << "Using an NPTF luminosity function" << std::endl << std::endl;
        nptfLBreak = NPTF_L_BREAK;
        nptfLMin = 0;
        break;
    case LUMINOSITY_FUNCTION::BARTELS_15:
        std::cout << "Using the wavelet fit luminosity function" << std::endl << std::endl;
        break;
    case LUMINOSITY_FUNCTION::BARTELS_18:
        std::cout << "Using the Bartels 18 luminosity function" << std::endl << std::endl;
        nptfLBreak = BARTELS_18_L_BREAK;
        nptfLMin = 1e29;
        break;
    }


    for (double i = 0; i < FLUX_HIST_BINS; i++) {
        double low = pow(10, (log10(FLUX_HIST_HIGH) - log10(FLUX_HIST_LOW)) * (i / FLUX_HIST_BINS) + log10(FLUX_HIST_LOW));
        double high = pow(10, (log10(FLUX_HIST_HIGH) - log10(FLUX_HIST_LOW)) * ((i + 1) / FLUX_HIST_BINS) + log10(FLUX_HIST_LOW));
        fluxValues[i] = (low + high) / 2;
        fluxWidths[i] = high - low;
        fluxLow[i] = low;
        fluxHigh[i] = high;
    }

    std::array<double, FLUX_HIST_BINS> probs;
    probs.fill(0);

    for (int x = 0; x < thresholds.skyShape.first; x++) {
        //std::cout << x << '/' << skyMap->size() << std::endl;
        for (int y = 0; y < thresholds.skyShape.second; y++) {
            CoordPair latLon = thresholds.indexToLatLon({ x, y });
            if (myAbs(latLon.first) > 2 * pi / 180) {// Cut out 2 degrees around the equator on each side.
                std::array<double, FLUX_HIST_BINS> probs_here = getHistAtLatLon(latLon, thresholds.getFluxThreshold(latLon));
                for (int i = 0; i < FLUX_HIST_BINS; i++) {
                    probs[i] += probs_here[i];
                }
            }
        }
    }

    double scale;
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, fermilabPaperPoint[1], fermilabPaperPoint[0]);
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalPaperPoint[0], logNormalPaperPoint[1]);
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalPloegPoint[0], logNormalPloegPoint[1]);
        break;
    case LUMINOSITY_FUNCTION::GAUTAM:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, logNormalGautamPoint[0], logNormalGautamPoint[1]);
        break;
    case LUMINOSITY_FUNCTION::NPTF:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, NPTF_N_1_HIST, NPTF_N_2_HIST);
        break;
    case LUMINOSITY_FUNCTION::BARTELS_15:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, BARTELS_15_L_MIN, bartels15Point[1]);
        break;
    case LUMINOSITY_FUNCTION::BARTELS_18:
        scale = FLUX_EXCESS / getValueAtConfig(VALUE::TOTAL_FLUX, BARTELS_18_N_1_HIST, BARTELS_18_N_2_HIST);
        break;
    };

    for (int i = 0; i < FLUX_HIST_BINS; i++) {
        probs[i] *= scale * CM_PER_KPC_SQUARED;//FLUX_EXCESS / probs[FLUX_HIST_BINS];
    }


    double sum = 0;
    for (int i = 0; i < FLUX_HIST_BINS; i++) {
        sum += probs[i];
    }
    std::cout << sum << std::endl;
    // Write data to an output file
    std::ofstream histFile;
#ifndef REAL_HISTOGRAM
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/power-law/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "alpha = " << ALPHA_HIST << ", L_min = " << fermilabPaperPoint[1] << ", L_max = " << fermilabPaperPoint[0] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/log-normal/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L0 = " << logNormalPaperPoint[0] << ", sigma = " << logNormalPaperPoint[1] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/ploeg/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L0 = " << logNormalPloegPoint[0] << ", sigma = " << logNormalPloegPoint[1] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::GAUTAM:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/gautam/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L0 = " << logNormalGautamPoint[0] << ", sigma = " << logNormalGautamPoint[1] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::NPTF:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/nptf/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L_break = " << NPTF_L_BREAK << ", n1 = " << NPTF_N_1_HIST << ", n2 = " << NPTF_N_2_HIST << std::endl;
        break;
    case LUMINOSITY_FUNCTION::BARTELS_15:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/wavelet-fit/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "header" << std::endl;
        break;
    case LUMINOSITY_FUNCTION::BARTELS_18:
        histFile.open(ROOT "luminosity-models-position/data-" + std::to_string((int)SENSITIVITY_DIVISOR) + "x/bartels18/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "header" << std::endl;
        break;
    }
#else
    switch (luminosityFunction) {
    case LUMINOSITY_FUNCTION::POWER_LAW:
    case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
        histFile.open(ROOT "luminosity-models-position/real-data/power-law/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "alpha = " << ALPHA_HIST << ", L_min = " << fermilabPaperPoint[1] << ", L_max = " << fermilabPaperPoint[0] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::LOG_NORMAL:
        histFile.open(ROOT "luminosity-models-position/real-data/log-normal/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L0 = " << logNormalPaperPoint[0] << ", sigma = " << logNormalPaperPoint[1] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::PLOEG:
        histFile.open(ROOT "luminosity-models-position/real-data/ploeg/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L0 = " << logNormalPloegPoint[0] << ", sigma = " << logNormalPloegPoint[1] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::GAUTAM:
        histFile.open(ROOT "luminosity-models-position/real-data/gautam/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L0 = " << logNormalGautamPoint[0] << ", sigma = " << logNormalGautamPoint[1] << std::endl;
        break;
    case LUMINOSITY_FUNCTION::NPTF:
        histFile.open(ROOT "luminosity-models-position/real-data/nptf/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "L_break = " << NPTF_L_BREAK << ", n1 = " << NPTF_N_1_HIST << ", n2 = " << NPTF_N_2_HIST << std::endl;
        break;
    case LUMINOSITY_FUNCTION::BARTELS_15:
        histFile.open(ROOT "luminosity-models-position/real-data/wavelet-fit/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "header" << std::endl;
        break;
    case LUMINOSITY_FUNCTION::BARTELS_18:
        histFile.open(ROOT "luminosity-models-position/real-data/bartels18/flux-hist-DR"
        + std::to_string((int)DATA_RELEASE) + ".txt");
        histFile << "header" << std::endl;
        break;
    }
#endif
    for (int i = 0; i < FLUX_HIST_BINS; i++) {
        histFile << probs[i] << ',';
    }
    histFile << std::endl;

    for (double value : fluxValues) {
        histFile << value << ',';
    }
    histFile << std::endl;

    std::cout << "Done." << std::endl;
    return 1;
}


int main(int argc, char** argv) {
    std::cout << "Sensitivity " << SENSITIVITY_DIVISOR << ", DR " << DATA_RELEASE << std::endl;
    METHOD method = METHOD::HIST;
    if (argc == 3) {
        if (strcmp(argv[2], "powerlaw") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::POWER_LAW;
        }
        else if (strcmp(argv[2], "lognormal") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::LOG_NORMAL;
        }
        else if (strcmp(argv[2], "ploeg") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::PLOEG;
        }
        else if (strcmp(argv[2], "gautam") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::GAUTAM;
        }
        else if (strcmp(argv[2], "nptf") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::NPTF;
        }
        else if (strcmp(argv[2], "alpha") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::POWER_LAW_ALPHA;
        }
        else if (strcmp(argv[2], "wavelet") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::BARTELS_15;
        }
        else if (strcmp(argv[2], "18") == 0) {
            luminosityFunction = LUMINOSITY_FUNCTION::BARTELS_18;
        }
        else {
            luminosityFunction = LUMINOSITY_FUNCTION::ERROR;
        }

        if (strcmp(argv[1], "count") == 0) {
            method = METHOD::COUNT;
        }
        else if (strcmp(argv[1], "hist") == 0) {
            method = METHOD::HIST;
        }
        else {
            method = METHOD::ERROR;
        }
    }

    if (method == METHOD::COUNT) {
        std::cout << "Count" << std::endl;
        switch (luminosityFunction) {
        case LUMINOSITY_FUNCTION::POWER_LAW:
            return powerLaw();
        case LUMINOSITY_FUNCTION::LOG_NORMAL:
            return logNormal();
        case LUMINOSITY_FUNCTION::PLOEG:
            return ploeg();
        case LUMINOSITY_FUNCTION::NPTF:
            return nptf();
        case LUMINOSITY_FUNCTION::POWER_LAW_ALPHA:
            return powerLawAlpha();
        case LUMINOSITY_FUNCTION::ERROR:
            std::cout << "The argument \"" + std::string(argv[2]) + "\" is not supported. Options are \"powerlaw\", \"lognormal\", \"ploeg\", \"alpha\" and \"nptf\"." << std::endl;
            return 0;
        }
    }

    else if (method == METHOD::HIST) {
        std::cout << "Hist" << std::endl;
        switch (luminosityFunction) {
        case LUMINOSITY_FUNCTION::ERROR:
            std::cout << "The argument \"" + std::string(argv[2]) + "\" is not supported. Options are \"powerlaw\", \"lognormal\", \"ploeg\", \"alpha\" and \"nptf\"." << std::endl;
            return 0;
        default:
            return hist();
        }
    }

    else {
        std::cout << "The argument \"" + std::string(argv[1]) + "\" is not supported. Options are \"count\" and \"hist\"." << std::endl;
        return 0;
    }
    return 1;
}
