
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

//ref:  https://github.com/mackron/dr_libs/blob/master/dr_wav.h
#define DR_WAV_IMPLEMENTATION

#include "dr_wav.h"
#include "PitchShift.h"
#include "timing.h"

void wavWrite_s16(char *filename, int16_t *buffer, int sampleRate, uint32_t totalSampleCount) {
    drwav_data_format format;
    format.container = drwav_container_riff;     // <-- drwav_container_riff = normal WAV files, drwav_container_w64 = Sony Wave64.
    format.format = DR_WAVE_FORMAT_PCM;          // <-- Any of the DR_WAVE_FORMAT_* codes.
    format.channels = 1;
    format.sampleRate = (drwav_uint32) sampleRate;
    format.bitsPerSample = 16;
    drwav *pWav = drwav_open_file_write(filename, &format);
    if (pWav) {
        drwav_uint64 samplesWritten = drwav_write(pWav, totalSampleCount, buffer);
        drwav_uninit(pWav);
        if (samplesWritten != totalSampleCount) {
            fprintf(stderr, "ERROR\n");
            exit(1);
        }
    }
}

int16_t *wavRead_s16(char *filename, uint32_t *sampleRate, uint64_t *totalSampleCount, uint32_t *channels) {
    int16_t *buffer = drwav_open_and_read_file_s16(filename, channels, sampleRate, totalSampleCount);
    if (buffer == NULL) {
        printf("ERROR.");
    }
    if (*channels == 2 && buffer != NULL) {
        int16_t *bufferSave = buffer;
        for (uint64_t i = 0; i < *totalSampleCount; i += 2) {
            *bufferSave++ = ((buffer[i] + buffer[i + 1]) >> 1);
        }
        *totalSampleCount = *totalSampleCount >> 1;
    } else if (*channels != 1) {
        drwav_free(buffer);
        buffer = NULL;
        *sampleRate = 0;
        *totalSampleCount = 0;
    }
    return buffer;
}

void splitpath(const char *path, char *drv, char *dir, char *name, char *ext) {
    const char *end;
    const char *p;
    const char *s;
    if (path[0] && path[1] == ':') {
        if (drv) {
            *drv++ = *path++;
            *drv++ = *path++;
            *drv = '\0';
        }
    } else if (drv)
        *drv = '\0';
    for (end = path; *end && *end != ':';)
        end++;
    for (p = end; p > path && *--p != '\\' && *p != '/';)
        if (*p == '.') {
            end = p;
            break;
        }
    if (ext)
        for (s = end; (*ext = *s++);)
            ext++;
    for (p = end; p > path;)
        if (*--p == '\\' || *p == '/') {
            p++;
            break;
        }
    if (name) {
        for (s = p; s < end;)
            *name++ = *s++;
        *name = '\0';
    }
    if (dir) {
        for (s = path; s < p;)
            *dir++ = *s++;
        *dir = '\0';
    }
}

int main(int argc, char *argv[]) {
    printf("Audio Processing \n");
    printf("blog:http://cpuimage.cnblogs.com/ \n");
    printf("Pitch Shifting Using The Fourier Transform\n");

    if (argc < 2)
        return -1;

    char *in_file = argv[1];
    uint32_t sampleRate = 0;
    uint64_t totalSampleCount = 0;
    uint32_t channels = 0;
    short *data_in = wavRead_s16(in_file, &sampleRate, &totalSampleCount, &channels);
    if (data_in != NULL) {
        float pitchShift = 0.9f;
        size_t ms = 50;
        size_t frameSize = sampleRate * ms / 1000;
        frameSize += frameSize % 2;
        planData pitchPlanData = {0};
        double startTime = now();
        makePlanData(frameSize, sampleRate, &pitchPlanData);
        pitchshift(pitchShift, data_in, data_in, totalSampleCount, &pitchPlanData);
        // turn to minion pitch
        {
            totalSampleCount /= 2;
            short *samples = data_in;
            for (int i = 0; i < totalSampleCount; i++) {
                data_in[i] = samples[0];
                samples += 2;
            }
        }
        double time_interval = calcElapsed(startTime, now());
        freePlanData(&pitchPlanData);
        printf("time interval: %f ms\n ", (time_interval * 1000));
    }
    char drive[3];
    char dir[256];
    char fname[256];
    char ext[256];
    char out_file[1024];
    splitpath(in_file, drive, dir, fname, ext);
    sprintf(out_file, "%s%s%s_out%s", drive, dir, fname, ext);
    wavWrite_s16(out_file, data_in, sampleRate, totalSampleCount);
    if (data_in) {
        free(data_in);
    }
    printf("press any key to exit.\n");
    getchar();
    return 0;
}
