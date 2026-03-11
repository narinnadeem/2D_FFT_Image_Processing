#define PI 3.14159265358979323846f
#include "math.h"

float mysin(float x);
float mycos(float x);
float myexp(float x);
int log2_int(int n);

float mysin(float x) {
    // 1. Wrap x to range [-PI, PI] for accuracy
    while (x > PI) x -= 2 * PI;
    while (x < -PI) x += 2 * PI;

    float res = 0.0f;
    float term = x; 
    float x_sq = x * x;
    
    for (int i = 1; i <= 10; i++) {
        res += term;
        // Calculate next term: multiply by -x^2 and divide by (2i)(2i+1)
        term *= -x_sq / ((2 * i) * (2 * i + 1));
    }
    return res;
}

float mycos(float x) {
    return mysin(x + (PI / 2.0f));
}

float myexp(float x) {
    if (x < 0) return 1.0f / myexp(-x); // Handle negative exponents

    float res = 1.0f;
    float term = 1.0f;
    for (int i = 1; i <= 25; i++) { 
        term *= x / (float)i;
        res += term;
    }
    return res;
}

int log2_int(int n) {
    int log = 0;
    while (n > 1) {
        n >>= 1; 
        log++;
    }
    return log;
}
