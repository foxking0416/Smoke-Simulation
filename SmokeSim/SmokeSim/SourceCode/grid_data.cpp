#include "grid_data.h"


GridData::GridData() :
   mDfltValue(0.0), mMax(0.0,0.0,0.0)
{
}

GridData::GridData(const GridData& orig) :
   mDfltValue(orig.mDfltValue)
{
   mData = orig.mData;
   mMax = orig.mMax;
}

GridData::~GridData() 
{
}

std::vector<double>& GridData::data()
{
   return mData;
}

GridData& GridData::operator=(const GridData& orig)
{
   if (this == &orig)
   {
      return *this;
   }
   mDfltValue = orig.mDfltValue;
   mData = orig.mData;
   mMax = orig.mMax;
   return *this;
}

void GridData::initialize(double dfltValue)
{
   mDfltValue = dfltValue;
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridData::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

const double GridData::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

void GridData::getCell(const vec3& pt, int& i, int& j, int& k)
{
   vec3 pos = worldToSelf(pt); 
   i = (int) (pos[0]/theCellSize);
   j = (int) (pos[1]/theCellSize);
   k = (int) (pos[2]/theCellSize);   
}

double GridData::CERT(double valueM1, double value, double valueP1, double valueP2, double t)
{

	double deltaK = valueP1 - value;
	double dK = (valueP1 - valueM1) / 2.0f;
	double dKPlus1 = (valueP2 - value) / 2.0f;
	
	if(deltaK == 0)
	{
		dK = 0;
		dKPlus1 = 0;
	}
	else if(deltaK * dK < 0)
	{
		dK = 0;
	}
	else if(deltaK * dKPlus1 < 0)
	{
		dKPlus1 = 0;
	}


	double a0 = value;
	double a1 = dK;
	double a2 = 3*deltaK - 2*dK - dKPlus1;
	double a3 = dK + dKPlus1 - 2*deltaK;

	double f = a3 * pow(t,3) + a2 * pow(t,2) + a1 * t + a0;

	return f;
}
	
double GridData::interpolate(const vec3& pt)
{

	// TODO: Implement sharper cubic interpolation here.

	vec3 pos = worldToSelf(pt);

	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);
	int k = (int) (pos[2]/theCellSize);

	double scale = 1.0/theCellSize;  
	double fractx = scale*(pos[0] - i*theCellSize);
	double fracty = scale*(pos[1] - j*theCellSize);
	double fractz = scale*(pos[2] - k*theCellSize);

	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);

	//// LINEAR INTERPOLATION:
	//	// Y @ low X, low Z:
	//double tmp1 = (*this)(i,j,k);
	//double tmp2 = (*this)(i,j+1,k);
	//// Y @ high X, low Z:
	//double tmp3 = (*this)(i+1,j,k);
	//double tmp4 = (*this)(i+1,j+1,k);

	//// Y @ low X, high Z:
	//double tmp5 = (*this)(i,j,k+1);
	//double tmp6 = (*this)(i,j+1,k+1);
	//// Y @ high X, high Z:
	//double tmp7 = (*this)(i+1,j,k+1);
	//double tmp8 = (*this)(i+1,j+1,k+1);

	//// Y @ low X, low Z
	//double tmp12 = LERP(tmp1, tmp2, fracty);
	//// Y @ high X, low Z
	//double tmp34 = LERP(tmp3, tmp4, fracty);

	//// Y @ low X, high Z
	//double tmp56 = LERP(tmp5, tmp6, fracty);
	//// Y @ high X, high Z
	//double tmp78 = LERP(tmp7, tmp8, fracty);

	//// X @ low Z
	//double tmp1234 = LERP (tmp12, tmp34, fractx);
	//// X @ high Z
	//double tmp5678 = LERP (tmp56, tmp78, fractx);

	//// Z
	//double tmp = LERP(tmp1234, tmp5678, fractz);



	double atmp1  = (*this)(i-1, j-1, k-1);
	double atmp2  = (*this)(i  , j-1, k-1);
	double atmp3  = (*this)(i+1, j-1, k-1);
	double atmp4  = (*this)(i+2, j-1, k-1);
	double atmp5  = (*this)(i-1, j  , k-1);
	double atmp6  = (*this)(i  , j  , k-1);
	double atmp7  = (*this)(i+1, j  , k-1);
	double atmp8  = (*this)(i+2, j  , k-1);
	double atmp9  = (*this)(i-1, j+1, k-1);
	double atmp10 = (*this)(i  , j+1, k-1);
	double atmp11 = (*this)(i+1, j+1, k-1);
	double atmp12 = (*this)(i+2, j+1, k-1);
	double atmp13 = (*this)(i-1, j+2, k-1);
	double atmp14 = (*this)(i  , j+2, k-1);
	double atmp15 = (*this)(i+1, j+2, k-1);
	double atmp16 = (*this)(i+2, j+2, k-1);
									  
	double atmp17 = (*this)(i-1, j-1, k);
	double atmp18 = (*this)(i  , j-1, k);
	double atmp19 = (*this)(i+1, j-1, k);
	double atmp20 = (*this)(i+2, j-1, k);
	double atmp21 = (*this)(i-1, j  , k);
	double atmp22 = (*this)(i  , j  , k);
	double atmp23 = (*this)(i+1, j  , k);
	double atmp24 = (*this)(i+2, j  , k);
	double atmp25 = (*this)(i-1, j+1, k);
	double atmp26 = (*this)(i  , j+1, k);
	double atmp27 = (*this)(i+1, j+1, k);
	double atmp28 = (*this)(i+2, j+1, k);
	double atmp29 = (*this)(i-1, j+2, k);
	double atmp30 = (*this)(i  , j+2, k);
	double atmp31 = (*this)(i+1, j+2, k);
	double atmp32 = (*this)(i+2, j+2, k);

	double atmp33 = (*this)(i-1, j-1, k+1);
	double atmp34 = (*this)(i  , j-1, k+1);
	double atmp35 = (*this)(i+1, j-1, k+1);
	double atmp36 = (*this)(i+2, j-1, k+1);
	double atmp37 = (*this)(i-1, j  , k+1);
	double atmp38 = (*this)(i  , j  , k+1);
	double atmp39 = (*this)(i+1, j  , k+1);
	double atmp40 = (*this)(i+2, j  , k+1);
	double atmp41 = (*this)(i-1, j+1, k+1);
	double atmp42 = (*this)(i  , j+1, k+1);
	double atmp43 = (*this)(i+1, j+1, k+1);
	double atmp44 = (*this)(i+2, j+1, k+1);
	double atmp45 = (*this)(i-1, j+2, k+1);
	double atmp46 = (*this)(i  , j+2, k+1);			  
	double atmp47 = (*this)(i+1, j+2, k+1);
	double atmp48 = (*this)(i+2, j+2, k+1);

	double atmp49 = (*this)(i-1, j-1, k+2);
	double atmp50 = (*this)(i  , j-1, k+2);
	double atmp51 = (*this)(i+1, j-1, k+2);
	double atmp52 = (*this)(i+2, j-1, k+2);
	double atmp53 = (*this)(i-1, j  , k+2);
	double atmp54 = (*this)(i  , j  , k+2);
	double atmp55 = (*this)(i+1, j  , k+2);
	double atmp56 = (*this)(i+2, j  , k+2);
	double atmp57 = (*this)(i-1, j+1, k+2);
	double atmp58 = (*this)(i  , j+1, k+2);
	double atmp59 = (*this)(i+1, j+1, k+2);
	double atmp60 = (*this)(i+2, j+1, k+2);
	double atmp61 = (*this)(i-1, j+2, k+2);
	double atmp62 = (*this)(i  , j+2, k+2);
	double atmp63 = (*this)(i+1, j+2, k+2);
	double atmp64 = (*this)(i+2, j+2, k+2);


	double temp_1_5_9_13_Y1  = CERT(atmp1, atmp5, atmp9,  atmp13, fracty);
	double temp_2_6_10_14_Y2 = CERT(atmp2, atmp6, atmp10, atmp14, fracty);
	double temp_3_7_11_15_Y3 = CERT(atmp3, atmp7, atmp11, atmp15, fracty);
	double temp_4_8_12_16_Y4 = CERT(atmp4, atmp8, atmp12, atmp16, fracty);

	double temp_17_21_25_29_Y5 = CERT(atmp17, atmp21, atmp25, atmp29, fracty);
	double temp_18_22_26_30_Y6 = CERT(atmp18, atmp22, atmp26, atmp30, fracty);
	double temp_19_23_27_31_Y7 = CERT(atmp19, atmp23, atmp27, atmp31, fracty);
	double temp_20_24_28_32_Y8 = CERT(atmp20, atmp24, atmp28, atmp32, fracty);

	double temp_33_37_41_45_Y9  = CERT(atmp33, atmp37, atmp41, atmp45, fracty);
	double temp_34_38_42_46_Y10 = CERT(atmp34, atmp38, atmp42, atmp46, fracty);
	double temp_35_39_43_47_Y11 = CERT(atmp35, atmp39, atmp43, atmp47, fracty);
	double temp_36_40_44_48_Y12 = CERT(atmp36, atmp40, atmp44, atmp48, fracty);

	double temp_49_53_57_61_Y13 = CERT(atmp49, atmp53, atmp57, atmp61, fracty);
	double temp_50_54_58_62_Y14 = CERT(atmp50, atmp54, atmp58, atmp62, fracty);
	double temp_51_55_59_63_Y15 = CERT(atmp51, atmp55, atmp59, atmp63, fracty);
	double temp_52_56_60_64_Y16 = CERT(atmp52, atmp56, atmp60, atmp64, fracty);

	double temp_Y1_Y2_Y3_Y4 = CERT(temp_1_5_9_13_Y1,     temp_2_6_10_14_Y2,    temp_3_7_11_15_Y3,    temp_4_8_12_16_Y4,    fractx);
	double temp_Y5_Y6_Y7_Y8 = CERT(temp_17_21_25_29_Y5,  temp_18_22_26_30_Y6,  temp_19_23_27_31_Y7,  temp_20_24_28_32_Y8,  fractx);
	double temp_Y9_Y10_Y11_Y12 = CERT(temp_33_37_41_45_Y9,  temp_34_38_42_46_Y10, temp_35_39_43_47_Y11, temp_36_40_44_48_Y12, fractx);
	double temp_Y13_Y14_Y15_Y16 = CERT(temp_49_53_57_61_Y13, temp_50_54_58_62_Y14, temp_51_55_59_63_Y15, temp_52_56_60_64_Y16, fractx);

	double atmp =  CERT(temp_Y1_Y2_Y3_Y4, temp_Y5_Y6_Y7_Y8, temp_Y9_Y10_Y11_Y12, temp_Y13_Y14_Y15_Y16, fractz);

	return atmp;
	
}

vec3 GridData::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
   return out;
}

GridDataX::GridDataX() : GridData()
{
}

GridDataX::~GridDataX()
{
}

void GridDataX::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*(theDim[0]+1);
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize((theDim[0]+1)*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataX::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

const double GridDataX::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataX::worldToSelf(const vec3& pt) const
{   
   vec3 out;
   out[0] = min(max(0.0, pt[0]), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataY::GridDataY() : GridData()
{
}

GridDataY::~GridDataY()
{
}

void GridDataY::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*(theDim[1]+1);
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*(theDim[1]+1)*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataY::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

const double GridDataY::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataY::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataZ::GridDataZ() : GridData()
{
}

GridDataZ::~GridDataZ()
{
}

void GridDataZ::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*(theDim[2]+1);
   mData.resize(theDim[0]*theDim[1]*(theDim[2]+1), false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataZ::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

const double GridDataZ::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

vec3 GridDataZ::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]), mMax[2]);
   return out;
}
