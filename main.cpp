#include "processor.hpp"

// ============================================================================

// #define _SVGMRF			// ���ꕪ�U���琄�肷��摜�f�[�^GMRF
// #define _UpScale		// GMRF��p��������
// #define _DVGMRF			// �������U���琄�肷��摜�f�[�^GMRF
// #define _SVODGMRF		// 1�����f�[�^�ɑ΂���GMRF
// #define _DVODGMRF		// 1�����f�[�^�ɑ΂���DVGMRF

// #define _DVSVCOMP		// DVGMRF��SVGMRF�̔�r����(�摜)
// #define _DVSVODCOMP		// DVODGMRF��SVODGMRF�̔�r����(1�����f�[�^)

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