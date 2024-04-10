#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void c_set_seedname(char c_seedname[80]);
void c_read_ddp();
double c_ddp_rcut();
void c_eval_ddp(
    const int num_ions,
    char* ion_names,
    const double *ion_positions,
    const double *lattice_car,
    double *e,
    double *f,
    double *s);

#ifdef __cplusplus
}
#endif
