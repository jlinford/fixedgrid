#ifndef __SAPRC99_H__
#define __SAPRC99_H__

void init_saprc99();
void init_conc(double *conc);
void saprc99(double *conc, double temp, double tstart, double tend, double dt);

#endif
