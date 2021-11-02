#include <limits.h>
#include "fix16.h"

fix16_t fix16_sin_parabola(fix16_t inAngle)
{
	fix16_t abs_inAngle, abs_retval, retval;
	fix16_t mask;

	/* Absolute function */
	mask = (inAngle >> (sizeof(fix16_t)*CHAR_BIT-1));
	abs_inAngle = (inAngle + mask) ^ mask;
	
	/* On 0->PI, sin looks like x² that is :
	   - centered on PI/2,
	   - equals 1 on PI/2,
	   - equals 0 on 0 and PI
	  that means :  4/PI * x  - 4/PI² * x²
	  Use abs(x) to handle (-PI) -> 0 zone.
	 */
	retval = fix16_mul(FOUR_DIV_PI, inAngle) + fix16_mul( fix16_mul(_FOUR_DIV_PI2, inAngle), abs_inAngle );
	/* At this point, retval equals sin(inAngle) on important points ( -PI, -PI/2, 0, PI/2, PI),
	   but is not very precise between these points
	 */

	return retval;
}

fix16_t fix16_sin(fix16_t inAngle) {
	fix16_t tempAngle = inAngle % (fix16_pi << 1);

	
	if(tempAngle > fix16_pi)
		tempAngle -= (fix16_pi << 1);
	else if(tempAngle < -fix16_pi)
		tempAngle += (fix16_pi << 1);

	

	fix16_t tempAngleSq = fix16_mul(tempAngle, tempAngle);

	
	fix16_t tempOut;
	tempOut = fix16_mul(-13, tempAngleSq) + 546;
	tempOut = fix16_mul(tempOut, tempAngleSq) - 10923;
	tempOut = fix16_mul(tempOut, tempAngleSq) + 65536;
	tempOut = fix16_mul(tempOut, tempAngle);




	return tempOut;
}

fix16_t fix16_cos(fix16_t inAngle) {
	return fix16_sin(inAngle + (fix16_pi >> 1));
}

fix16_t fix16_tan(fix16_t inAngle) {
	return fix16_sdiv(fix16_sin(inAngle), fix16_cos(inAngle));
}

fix16_t fix16_asin(fix16_t inValue) {
	if((inValue > fix16_one) || (inValue < -fix16_one))
		return 0;
	fix16_t tempOut;
	tempOut = (fix16_one - fix16_mul(inValue, inValue));
	tempOut = fix16_div(inValue, fix16_sqrt(tempOut));
	tempOut = fix16_atan(tempOut);
	return tempOut;
}

fix16_t fix16_acos(fix16_t inValue) {
	return ((fix16_pi >> 1) - fix16_asin(inValue));
}

fix16_t fix16_atan2(fix16_t inY , fix16_t inX) {
	fix16_t abs_inY, mask, angle, r, r_3;



	/* Absolute inY */
	mask = (inY >> (sizeof(fix16_t)*CHAR_BIT-1));
	abs_inY = (inY + mask) ^ mask;

	if (inX >= 0)
	{
		r = fix16_div( (inX - abs_inY), (inX + abs_inY));
		r_3 = fix16_mul(fix16_mul(r, r),r);
		angle = fix16_mul(0x00003240 , r_3) - fix16_mul(0x0000FB50,r) + PI_DIV_4;
	} else {
		r = fix16_div( (inX + abs_inY), (abs_inY - inX));
		r_3 = fix16_mul(fix16_mul(r, r),r);
		angle = fix16_mul(0x00003240 , r_3) - fix16_mul(0x0000FB50,r) + THREE_PI_DIV_4;
	}
	if (inY < 0)
	{
		angle = -angle;
	}


	return angle;
}

fix16_t fix16_atan(fix16_t inValue) {
	return fix16_atan2(inValue, fix16_one);
}
