// Copyright 2017 Dirk Hornung

#include "./constants.h"

double Constants::pifac = 24.*pow(M_PI*Vud*fpi, 2)*SEW;
double Constants::dpifac = pifac*sqrt(4.*pow(dVud/Vud, 2)
                                      +pow(dSEW/SEW, 2)+4.*pow(dfpi/fpi, 2));




