#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#ifdef USE_DOUBLE_TYPE
#define sample_real_t double
#else
#define sample_real_t  float
#endif

#define STB_FFT_IMPLEMENTAION

#include "stb_fft.h"

typedef struct {
    sample_real_t *inFifo;
    sample_real_t *outFifo;
    sample_real_t *fftWorkspace;
    sample_real_t *lastPhase;
    sample_real_t *phaseSum;
    sample_real_t *windowing;
    sample_real_t *outputAccumulator;
    sample_real_t *synthesizedFrequency;
    sample_real_t *synthesizedMagnitude;
    size_t frameSize;
    size_t overSampling;
    int64_t overlap;
    int64_t sampleRate;
} planData;


void freePlanData(planData *data) {
    if (data->inFifo)
        free(data->inFifo);
    data->inFifo = 0;
    if (data->outFifo)
        free(data->outFifo);
    data->outFifo = 0;
    if (data->fftWorkspace)
        free(data->fftWorkspace);
    data->fftWorkspace = 0;
    if (data->windowing)
        free(data->windowing);
    data->windowing = 0;
    if (data->lastPhase)
        free(data->lastPhase);
    data->lastPhase = 0;
    if (data->phaseSum)
        free(data->phaseSum);
    data->phaseSum = 0;
    if (data->outputAccumulator)
        free(data->outputAccumulator);
    data->outputAccumulator = 0;
    if (data->synthesizedFrequency)
        free(data->synthesizedFrequency);
    data->synthesizedFrequency = 0;
    if (data->synthesizedMagnitude)
        free(data->synthesizedMagnitude);
    data->synthesizedMagnitude = 0;
}

int makePlanData(size_t frameSize, size_t sampleRate, planData *data) {
    data->inFifo = (sample_real_t *) calloc(frameSize, sizeof(sample_real_t));
    data->outFifo = (sample_real_t *) calloc(frameSize, sizeof(sample_real_t));
    data->windowing = (sample_real_t *) calloc(frameSize, sizeof(sample_real_t));
    data->fftWorkspace = (sample_real_t *) calloc(frameSize * 2, sizeof(sample_real_t));
    data->lastPhase = (sample_real_t *) calloc(frameSize / 2 + 1, sizeof(sample_real_t));
    data->phaseSum = (sample_real_t *) calloc(frameSize / 2 + 1, sizeof(sample_real_t));
    data->outputAccumulator = (sample_real_t *) calloc(frameSize * 2, sizeof(sample_real_t));
    data->synthesizedFrequency = (sample_real_t *) calloc(frameSize, sizeof(sample_real_t));
    data->synthesizedMagnitude = (sample_real_t *) calloc(frameSize, sizeof(sample_real_t));
    if (!data->inFifo ||
        !data->outFifo ||
        !data->fftWorkspace ||
        !data->windowing ||
        !data->lastPhase ||
        !data->phaseSum ||
        !data->outputAccumulator ||
        !data->synthesizedFrequency ||
        !data->synthesizedMagnitude) {
        freePlanData(data);
        return 0;
    }
    data->frameSize = frameSize;
    data->sampleRate = sampleRate;
    return 1;
}

void pitchshift(sample_real_t pitch, const short *in, short *out, size_t samples, planData *data) {
    if (data == NULL)
        return;
    size_t frameSize = data->frameSize;
    size_t overSampling = 4;
    sample_real_t *inFifo = data->inFifo;
    sample_real_t *outFifo = data->outFifo;
    sample_real_t *fftWorkspace = data->fftWorkspace;
    sample_real_t *lastPhase = data->lastPhase;
    sample_real_t *phaseSum = data->phaseSum;
    sample_real_t *windowing = data->windowing;
    sample_real_t *outputAccumulator = data->outputAccumulator;
    sample_real_t *synthesizedFrequency = data->synthesizedFrequency;
    sample_real_t *synthesizedMagnitude = data->synthesizedMagnitude;
    sample_real_t magnitude, phase, deltaPhase, real, imag;
    int64_t i, k, qpd, index;
    int64_t halfFrameSize = frameSize / 2 + 1;
    int64_t step = frameSize / overSampling;
    sample_real_t binFrequencies = (sample_real_t) data->sampleRate / (sample_real_t) frameSize;
    sample_real_t expected = STB_TWOPI / overSampling;
    int64_t fifoLatency = frameSize - step;
    if (data->overlap == 0)
        data->overlap = fifoLatency;
    for (k = 0; k < frameSize; k++) {
        windowing[k] = -.5f * cosf(STB_TWOPI * k / frameSize) + 0.5f;
    }
    sample_real_t pitchWeight = pitch * binFrequencies;
    sample_real_t oversamp_weight = (overSampling / STB_TWOPI) * pitchWeight;
    sample_real_t meanExpected = expected / binFrequencies;
    stb_fft_real_plan *plan = NULL;
    int plan_bytes = stb_fft_real_plan_dft_1d(frameSize, NULL);
    if (plan_bytes > 0) {
        plan = (stb_fft_real_plan *) calloc(plan_bytes, 1);
        if (plan != NULL) {
            stb_fft_real_plan_dft_1d(frameSize, plan);
        }
    }
    if (plan) {
        for (i = 0; i < samples; i++) {
            inFifo[data->overlap] = in[i];
            out[i] = outFifo[data->overlap - fifoLatency];
            data->overlap++;
            if (data->overlap >= frameSize) {
                data->overlap = fifoLatency;
                for (k = 0; k < frameSize; k++) {
                    fftWorkspace[k] = inFifo[k] * windowing[k];
                }
                stb_fft_r2c_exec(plan, fftWorkspace, (cmplx *) fftWorkspace);
                cmplx *outCmplx = (cmplx *) fftWorkspace;
                memset(synthesizedMagnitude, 0, frameSize * sizeof(sample_real_t));
                memset(synthesizedFrequency, 0, frameSize * sizeof(sample_real_t));
                for (k = 0; k < halfFrameSize; k++) {
                    index = k * pitch;
                    if (index < halfFrameSize) {
                        real = outCmplx[k].real;
                        imag = outCmplx[k].imag;
                        magnitude = sqrtf(real * real + imag * imag);
                        phase = atan2f(imag, real);
                        deltaPhase = (phase - lastPhase[k]) - k * expected;
                        qpd = deltaPhase / STB_PI;
                        if (qpd >= 0)
                            qpd += qpd & 1;
                        else
                            qpd -= qpd & 1;
                        deltaPhase -= STB_PI * (sample_real_t) qpd;
                        lastPhase[k] = phase;
                        synthesizedMagnitude[index] += magnitude;
                        synthesizedFrequency[index] = k * pitchWeight + oversamp_weight * deltaPhase;
                    }
                }
                for (k = 0; k < halfFrameSize; k++) {
                    phaseSum[k] += meanExpected * synthesizedFrequency[k];
                    phase = phaseSum[k];
                    magnitude = synthesizedMagnitude[k];
                    stbSinCos(phase, &outCmplx[k].imag, &outCmplx[k].real);
                    outCmplx[k].real *= magnitude;
                    outCmplx[k].imag *= magnitude;
                }
                memset(fftWorkspace + (frameSize + 2), 0, sizeof(sample_real_t) * (frameSize - 2));
                stb_fft_c2r_exec(plan, (cmplx *) fftWorkspace, fftWorkspace);
                sample_real_t accOversamp = 2.f / (halfFrameSize * overSampling);
                for (k = 0; k < frameSize; k++) {
                    outputAccumulator[k] += windowing[k] * fftWorkspace[k] * accOversamp;
                }
                memcpy(outFifo, outputAccumulator, step * sizeof(sample_real_t));
                memmove(outputAccumulator, outputAccumulator + step, frameSize * sizeof(sample_real_t));
                memmove(inFifo, inFifo + step, fifoLatency * sizeof(sample_real_t));
            }
        }
        free(plan);
    }
}