#include "processor.hpp"

// ============================================================================

// #define _SVGMRF			// 同一分散から推定する画像データGMRF
// #define _UpScale		// GMRFを用いた超解像
// #define _DVGMRF			// 複数分散から推定する画像データGMRF
// #define _SVODGMRF		// 1次元データに対するGMRF
// #define _DVODGMRF		// 1次元データに対するDVGMRF

// #define _DVSVCOMP		// DVGMRFとSVGMRFの比較実験(画像)
// #define _DVSVODCOMP		// DVODGMRFとSVODGMRFの比較実験(1次元データ)

// #define _IVHGMRF
// #define _DVHGMRF
// #define _HGMRFCOMP

#define _GMRFCOMP

// ============================================================================

int main(int argc, char* argv[]) 
{
#ifdef _SVGMRF
	SVGMRF::processor svgmrf;
	svgmrf.process();
#endif // SVGMRF

#ifdef _UpScale
	GMRFUpScale::processor gmrfus;
	gmrfus.process();
#endif // UpScale

#ifdef _DVGMRF
	DVGMRF::processor dvgmrf;
	dvgmrf.process();
#endif // DVGMRF

#ifdef _SVODGMRF
	SVODGMRF::processor svodgmrf;
	svodgmrf.process();
#endif // SVODGMRF

#ifdef _DVODGMRF
	DVODGMRF::processor dvodgmrf;
	dvodgmrf.process();
#endif // DVODGMRF

#ifdef _DVSVCOMP
	DVSVCOMP::processor dvsvcomp;
	dvsvcomp.process();
#endif // DVSVCOMP

#ifdef _DVSVODCOMP
	DVSVODCOMP::processor dvsvodcomp;
	dvsvodcomp.process();
#endif // DVSVODCOMP

#ifdef _IVHGMRF
	IVHGMRF::processor ivhgmrf;
	ivhgmrf.process();
#endif // IVHGMRF

#ifdef _DVHGMRF
	DVHGMRF::processor dvhgmrf;
	dvhgmrf.process();
#endif // _DVHGMRF

#ifdef _HGMRFCOMP
	HGMRFCOMP::processor hgmrf_comp;
	hgmrf_comp.process();
#endif // HGMRFCOMP

#ifdef _GMRFCOMP
	GMRFCOMP::processor gmrf_comp;
	gmrf_comp.process();
#endif // GMRFCOMP
	return 0;
}