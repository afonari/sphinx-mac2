sed -e'
s/SxDiracVec</SxNewVec</g
s/SxDiracMat</SxNewVec</g
s/SxVector</SxNewVec</g
s/SxMatrix</SxNewVec</g
s/<Double>/<double>/g
s/<Complex8>/<SxComplex8>/g
s/<Complex16>/<SxComplex16>/g
s/<Int>/<int>/g
s/TReal8/double/g
s/TPrecPhi/double/g
s/TPrecCoeffG/PrecCoeffG/g
s/TPrecTauR/PrecTauR/g
s/TPrecFFTIDx/PrecFFTIdx/g
s/typename T::Type/T/g
s/#include <SxDirac.h>/#include <SxNewVec.h>/
s/VALIDATE_VECTOR/SX_VALIDATE_NEWVECTOR/
s/\.nRows *()/\.getNRows ()/g
s/\.nCols *()/\.getNCols ()/g
s/->nRows *()/->getNRows ()/g
s/->nCols *()/->getNCols ()/g
s/\.handle->auxData/.auxData/g
'  "$@"
