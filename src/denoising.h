#ifndef _DENOISING_H
#define _DENOISING_H

void andtv_filter(float *x, int N, float tau, float lambda); 

void aldtv_filter(float *x, int N, float tau, float lambda);

void diffusion_filter(float *x, int N, int T);

void moving_average_filter(float *x, int N, int L);

void gaussian_filter(float *x, int N, float *w, int L);

#endif
