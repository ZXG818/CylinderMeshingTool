#include "udf.h"
#include "mem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FUEL_START 1
#define POWER_TOTAL 4700.0 // 20cm rod, power is 4700 W
#define FUEL_TOTAL 7840
#define WATER_THREAD_ID 4
#define CLAD_THREAD_ID 5
#define FUEL_THREAD_ID 6

#define PI               3.1415926
#define ERR              1E-6
#define BOARD_X_UP       0.6
#define BOARD_X_DOWN    -0.6
#define BOARD_Y_UP       0.6
#define BOARD_Y_DOWN    -0.6
#define OUTER_X_UP       0.75
#define OUTER_X_DOWN    -0.75
#define OUTER_Y_UP       0.75
#define OUTER_Y_DOWN    -0.75
#define OUTER_Z_DOWN     0.0
#define OUTER_Z_UP       20.0

#define WATER_TYPE       1
#define FUEL_TYPE        2
#define CLAD_TYPE        3
static double residual = 0.01;

// MACRO FUNCTION
#define AD2DEG(theta) ((theta)*PI / 180.0)
#define CIRCLE_FUN_Y(R,x) (sqrt((R)*(R) - (x)*(x)))

// A Cylinder Mesh Method
typedef struct _POINT_
{
  int cell_id;
  int database;
  //int lw_database;
  double rho;
  double mcnp_volume; // the volume of the MCNP cell
  union
  {
    double power;
    double temperature;
  }u;
  double xyz[ND_ND];
}POINT, *PPOINT;


// global parameter
PPOINT g_coarse_fuel = NULL;
unsigned int g_index = 0;

// allocate the global fuel point array
void InitializeGlobal()
{
  g_coarse_fuel = (POINT *)malloc(sizeof(POINT)*FUEL_TOTAL);
  if(g_coarse_fuel == NULL)
  {
    Message("[*] InitializeGlobal allocate the memory error...\n");
    exit(0);
  }
  // set the g_coarse_fuel to zero memory
  memset(g_coarse_fuel, 0, sizeof(POINT)*FUEL_TOTAL);
}

void DestroyGlobal()
{
  // free the global allocated memory
  free(g_coarse_fuel);
  g_coarse_fuel = NULL;
}


/************************************************************************/
/*                              IDW METHOD                              */
/************************************************************************/

// Use the IDW method to make the interpolation from Fluent to MCNP.
void IDW_TO_MCNP(int thread_id, 
  PPOINT target, int N)
{
  Domain *d;
  Thread *t;
  cell_t c;
  double xyz[ND_ND];
  int i, j;
  double distance = 0.0;
  double temperature = 0.0;
  double rho = 0.0;
  double denomerator = 0.0;
  double numerator = 0.0;
  double x, y, z;
  double temp_x, temp_y, temp_z;

  d = Get_Domain(1);
  t = Lookup_Thread(d, thread_id);
  
  for (i=0; i<N; i++)
  {
    denomerator = 0.0;   // distance_total
    numerator = 0.0;
    temperature = 0.0;
    rho = 0.0;
    
    /* loop over the cells on the current thread */
    begin_c_loop(c, t)
    {
      /* caution please! the UDF's unit is M, and the MCNP's unit is CM, 
      *  so we should to multiply the 100 of the fuel plate location which in fluent */
      C_CENTROID(xyz, c, t);     // get the position
      temp_x = 100 * xyz[0];
      temp_y = 100 * xyz[1];
      temp_z = 100 * xyz[2];
      x = temp_x - target[i].xyz[0];
      y = temp_y - target[i].xyz[1];
      z = temp_z - target[i].xyz[2];
      distance = sqrt(x*x + y*y + z*z);
      if(distance <= ERR)
      {
        Message("[*] Here is an interpolation point very near of the target point...\n");
        target[i].u.temperature = C_T(c, t); 
        target[i].rho = C_R(c, t) / 1000.0;     // be cautious with here density of g/cc and kg/cube m
        break;
      }
      denomerator += pow(distance, -3.0);
      numerator    = pow(distance, -3.0);
      temperature += numerator * C_T(c, t);
      rho         += numerator * C_R(c, t);
    }end_c_loop(c, t)
    target[i].u.temperature = temperature / denomerator;
    target[i].rho = rho / denomerator / 1000.0;         // be cautious with here density of g/cc and kg/cube m;
  }
}


// Use the IDW method to make the interpolation from MCNP to Fluent.
void IDW_TO_FLUENT(int thread_id, 
  PPOINT point, int N)
{
  Domain *d;
  Thread *t;
  cell_t c;
  double xyz[ND_ND];
  int i, j;
  double distance = 0.0;
  double power = 0.0;
  double denomerator = 0.0;
  double numerator = 0.0;
  double x, y, z;
  double temp_x, temp_y, temp_z;
  
  d = Get_Domain(1);
  t = Lookup_Thread(d, thread_id);
  
  begin_c_loop(c, t)   // target point
  {
    denomerator = 0.0;
    numerator = 0.0;
    power = 0.0;
    C_CENTROID(xyz, c, t);
    /* caution please! the UDF's unit is M, and the MCNP's unit is CM, 
     *  so we should to multiply the 100 of the fuel plate location which in fluent */
    temp_x = 100 * xyz[0];
    temp_y = 100 * xyz[1];
    temp_z = 100 * xyz[2];
    /* loop over the cells of the current MCNP part */
    for(i=0; i<N; i++)    // interpolation points
    {
      x = temp_x - point[i].xyz[0];
      y = temp_y - point[i].xyz[1];
      z = temp_z - point[i].xyz[2];
      distance = sqrt(x*x + y*y + z*z);
      if(distance <= 1E-4)
      {
        Message("[*] Here is an interpolation point very near of the target point...\n");
        C_UDMI(c, t, 0) = point[i].u.power;
        break;
      }
      denomerator += pow(distance, -6.0);
      numerator    = pow(distance, -6.0);
      power       += numerator * point[i].u.power;
    }
    C_UDMI(c, t, 0) = power / denomerator;
  }end_c_loop(c, t)
}


/************************************************************************/
/*                              CORE CODE                               */
/************************************************************************/

// the cylinder part (fuel + cladding) grid generated function.
void McnpCylinderPartGrid(
  double R_IN, double R_OUT, int layer_cnt,
  double begin_z, double end_z, int z_cnt,
  int N,
  int kind,                 // kind is used to determined whether the cylinder or the annulus
  double *face_r,
  double *face_z,
  double up_tgx,            // up location of the tangent x
  double down_tgx,          // down location of the tangent x
  double up_tgy,            // up location of the tangent y
  double down_tgy,          // down location of the tangent y
  PPOINT array,
  double x_offset,          // the offset on X-axis of the cylinder circle location
  double y_offset)          // the offset on Y-axis of the cylinder circle location
{
  int i, j, k;
  int key = 0;
  double dr = (R_OUT - R_IN) / layer_cnt;
  double dz = (end_z - begin_z) / z_cnt;
  double d_theta = AD2DEG(360.0 / N);
  double r = 0.0;

  // r direction
  for (i = 0; i < layer_cnt; i++)
  {
    face_r[i] = R_IN + (i + 1) * dr;
  }

  // z direction
  for (i = 0; i < z_cnt + 1; i++)
  {
    face_z[i] = begin_z + i*dz;
  }

  key = 0;
  // begin to grid the cylinder
  for (i = 0; i < z_cnt; i++)              // z direction
  {
    for (j = 0; j < layer_cnt; j++)      // r direction
    {
      for (k = 0; k < N; k++)          // xy direction
      {
        // if the circle which contains the central point.
        if (j == 0 && kind == 1)
        {
          r = face_r[j] / 2.0;
          // get the current MCNP cell volume, unit: cubic cm
          array[key].mcnp_volume = dz*PI*face_r[j] * face_r[j] / (double)N;
        }
        // the circle which not contains the central point or the annulus
        else
        {
          if (j == 0)
          {
            r = (face_r[j] + R_IN) / 2.0;
            array[key].mcnp_volume = dz*((PI*face_r[j] * face_r[j]) - (PI*R_IN * R_IN)) / (double)N;
          }
          else
          {
            r = (face_r[j] + face_r[j - 1]) / 2.0;
            // get the current MCNP cell volume, unit: cubic cm
            array[key].mcnp_volume = dz*((PI*face_r[j] * face_r[j]) - (PI*face_r[j - 1] * face_r[j - 1])) / (double)N;
          }
        }
        array[key].xyz[0] = x_offset + r * cos(((2 * k + 1) * d_theta) / 2.0);   // x location
        array[key].xyz[1] = y_offset + r * sin(((2 * k + 1) * d_theta) / 2.0);   // y location
        array[key].xyz[2] = (face_z[i] + face_z[i + 1]) / 2.0;                   // z location
        key++;
      }
    }
  }
}


// the corner part coolant grid generated function.
void McnpCornerPartGrid(double R_OUT,
  double up_tgx,
  double *face_z,
  double begin_z, double end_z,
  int z_cnt, 
  double x_offset,
  double y_offset,
  int CornerXcnt,
  PPOINT array)
{
  int i, j, k;
  int key = 0; // key is the counter of the corner region points.
  double dx, dz;

  dx = up_tgx / (double)CornerXcnt;
  dz = (end_z - begin_z) / z_cnt;

  // discrete the Z direction
  for (i = 0; i < z_cnt + 1; i++)
  {
    face_z[i] = begin_z + i*dz;
  }

  // first is the up-right corner
  for (k = 0; k < z_cnt; k++)           // z direction
  {
    for (j = 0; j < CornerXcnt; j++)    // y direction
    {
      for (i = j; i < CornerXcnt; i++)  // x direction
      {
        array[key].xyz[0] = x_offset + (dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + (CIRCLE_FUN_Y(R_OUT, dx*j) + CIRCLE_FUN_Y(R_OUT, dx*(j + 1))) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  //second is the up-left corner
  for (k = 0; k < z_cnt; k++)
  {
    for (j = 0; j < CornerXcnt; j++)
    {
      for (i = j; i < CornerXcnt; i++)
      {
        array[key].xyz[0] = x_offset + (-1.0)*(dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + (CIRCLE_FUN_Y(R_OUT, dx*j) + CIRCLE_FUN_Y(R_OUT, dx*(j + 1))) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  // third is the bottom-left corner
  for (k = 0; k < z_cnt; k++)
  {
    for (j = 0; j < CornerXcnt; j++)
    {
      for (i = j; i < CornerXcnt; i++)
      {
        array[key].xyz[0] = x_offset + (-1.0)*(dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + (-1.0)*(CIRCLE_FUN_Y(R_OUT, dx*j) + CIRCLE_FUN_Y(R_OUT, dx*(j + 1))) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  // fourth is the bottom-right corner
  for (k = 0; k < z_cnt; k++)
  {
    for (j = 0; j < CornerXcnt; j++)
    {
      for (i = j; i < CornerXcnt; i++)
      {
        array[key].xyz[0] = x_offset + (dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + (-1.0)*(CIRCLE_FUN_Y(R_OUT, dx*j) + CIRCLE_FUN_Y(R_OUT, dx*(j + 1))) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
}

// the outer rectangular region grid generation
void OuterMainRectPartGrid(int a_cnt, int b_cnt,
  double *face_z,
  double begin_z, double end_z, int z_cnt,
  double x_offset, double y_offset,
  PPOINT array)
{
  int i, j, k;
  int key = 0;
  double da;
  double db;
  double dz;

  da = (OUTER_X_UP - BOARD_X_UP) / (double)a_cnt;
  db = (BOARD_X_UP - BOARD_X_DOWN) / (double)b_cnt;
  dz = (end_z - begin_z) / (double)z_cnt;

  // discrete the Z direction
  for (i = 0; i < z_cnt + 1; i++)
  {
    face_z[i] = begin_z + i*dz;
  }
  
  // first is the right rectangular region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < b_cnt; j++){
      for (i = 0; i < a_cnt; i++){
        array[key].xyz[0] = x_offset + BOARD_X_UP + (da*i + da*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + BOARD_Y_DOWN + (db*j + db*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  // second is the up rectangular region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < a_cnt; j++){
      for (i = 0; i < b_cnt; i++){
        array[key].xyz[0] = x_offset + BOARD_X_DOWN + (db*i + db*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + BOARD_Y_UP + (da*j + da*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  // third is the left rectangular region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < b_cnt; j++){
      for (i = 0; i < a_cnt; i++){
        array[key].xyz[0] = x_offset + OUTER_X_DOWN + (da*i + da*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + BOARD_Y_DOWN + (db*j + db*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  // fourth is the bottom rectangular region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < a_cnt; j++){
      for (i = 0; i < b_cnt; i++){
        array[key].xyz[0] = x_offset + BOARD_X_DOWN + (db*i + db*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + OUTER_Y_DOWN + (da*j + da*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
}

// the outer corner rectangular region grid generation (actually, it is a square.)
void OuterCornerRectPartGrid(int cnt, 
  double *face_z,
  double begin_z, double end_z, int z_cnt,
  double x_offset, double y_offset,
  PPOINT array)
{
  int i, j, k;
  int key = 0;
  double dx, dy;
  double dz;
  
  dx = dy = (OUTER_X_UP - BOARD_X_UP) / (double)cnt;
  dz = (end_z - begin_z) / (double)z_cnt;

  // discrete the Z direction
  for (i = 0; i < z_cnt + 1; i++)
    face_z[i] = begin_z + i*dz;

  // first is the up right corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        array[key].xyz[0] = x_offset + BOARD_X_UP + (dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + BOARD_Y_UP + (dy*j + dy*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  //second is the up left corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        array[key].xyz[0] = x_offset + OUTER_X_DOWN + (dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + BOARD_Y_UP + (dy*j + dy*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  // third is the down left corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        array[key].xyz[0] = x_offset + OUTER_X_DOWN + (dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + OUTER_Y_DOWN + (dy*j + dy*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
  // fourth is the down right corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        array[key].xyz[0] = x_offset + BOARD_X_UP + (dx*i + dx*(i + 1)) / 2.0;
        array[key].xyz[1] = y_offset + OUTER_Y_DOWN + (dy*j + dy*(j + 1)) / 2.0;
        array[key].xyz[2] = (face_z[k] + face_z[k + 1]) / 2.0;
        key++;
      }
    }
  }
}


// surface card generation
int McnpCylinderSurfaceCardInp(
  const char *filename,
  double *face_r, int layer_cnt,
  double *face_z, int z_cnt,
  int N,
  int kind,       // determine the cylinder(with central point), cylinder annulus, and cladding
  double up_tgx,
  double down_tgx,
  double up_tgy,
  double down_tgy,
  int surface_start,
  int *outer_cylinder,
  int *x_plate,       // x-axis plate
  int *y_plate,       // y-axis plate
  double x_offset,    // the offset on X-axis of the cylinder circle location
  double y_offset)    // the offset on Y-axis of the cylinder circle location   
{
  int i;
  int surface_number = 1;
  // line direction
  double beta_check = 0.0;
  double d_beta_check = 360.0 / N;
  double beta = 0.0;
  double d_beta = AD2DEG(d_beta_check);
  double K;
  double D = 0.0;
  FILE *fp = fopen(filename, "a");
  if (fp == NULL)
  {
    Message("[*] McnpSurfaceCardInp : open file error...\n");
    return -1;
  }

  surface_number = 1;

  // z direction
  for (i = 0; i < z_cnt + 1; i++)
  {
    // write the file
    fprintf(fp, "%d pz %.4f\n", surface_start + surface_number, face_z[i]);
    surface_number++;
  }

  // r direction
  for (i = 0; i < layer_cnt; i++)
  {
    // write the file
    // parallel with the Z-axis cylinder
    fprintf(fp, "%d c/z %.4f %.4f %.4f\n", surface_start + surface_number, x_offset, y_offset, face_r[i]);
    surface_number++;
  }
  *outer_cylinder = surface_start + surface_number - 1;

  // from X-axis to Y-axis  ---  anti-clockwise direction to make the division
  for (i = 0; i < N + 1; i++)
  {
    if (fabs(beta_check - 90.0) < ERR || fabs(beta_check - 270.0) < ERR)
    {
      D = x_offset;
      fprintf(fp, "%d p %.4f %.4f %.4f %.4f\n", surface_start + surface_number, 1.0, 0.0, 0.0, D);
      beta += d_beta;
      beta_check += d_beta_check;
      *y_plate = surface_start + surface_number;
      surface_number++;
      continue;
    }
    // get the K, which is the gradient
    K = tan(beta);
    // get the D, which is the intercept
    D = y_offset - K*x_offset;
    // write the file
    fprintf(fp, "%d p %.4f %.4f %.4f %.4f\n", surface_start + surface_number, -K, 1.0, 0.0, D);
    if (fabs(beta_check - 0.0) < ERR || fabs(beta_check - 180.0) < ERR)
    {
      *x_plate = surface_start + surface_number;
    }
    beta += d_beta;
    beta_check += d_beta_check;
    surface_number++;
  }

  // close the file
  fclose(fp);
  return surface_number;
}

// cell card
void McnpCylinderCellCardInp(const char *filename,
  int layer_cnt,
  int z_cnt,
  int N,
  int kind,
  int surface_start, int cell_start,
  int outer_cylinder,
  int outer_cylinder_backup,
  int x_plate,
  int y_plate,
  int surface_sum,
  PPOINT array)
{
  int i, j, k;
  int key = 1;
  double d_theta = 360.0 / N;
  double theta_sum = 0.0;
  int coeff1 = 1;
  int coeff2 = 1;
  FILE *fp = fopen(filename, "a");
  if (fp == NULL)
  {
    Message("[*] McnpCellPartCardInp : open file error...\n");
    return;
  }
  for (i = 1 + surface_start; i < 1 + surface_start + z_cnt; i++)                    // z direction
  {
    for (j = 2 + surface_start + z_cnt; j < 2 + surface_start + z_cnt + layer_cnt; j++)    // r direction
    {
      theta_sum = 0.0;
      for (k = 2 + surface_start + z_cnt + layer_cnt; k < 2 + surface_start + z_cnt + layer_cnt + N; k++) // xy direction
      {
        theta_sum += d_theta;
        if (theta_sum >= 0.0 && theta_sum <= 90.0 - d_theta)      { coeff1 =  1; coeff2 = -1; }
        else if (fabs(theta_sum - 90.0) <= ERR)                   { coeff1 =  1; coeff2 =  1; }
        else if (theta_sum > 90.0 && theta_sum <= 270.0 - d_theta){ coeff1 = -1; coeff2 =  1; }
        else if (fabs(theta_sum - 270.0) <= ERR)                  { coeff1 = -1; coeff2 = -1; }
        else if (theta_sum > 270.0 && theta_sum <= 360)           { coeff1 =  1; coeff2 = -1; }
        if (j == 2 + surface_start + z_cnt && kind == 1)  // the first inner cylinder radius direction
        {
          // cell_number material_number density R-DIR XY-DIR Z-DIR imp:n=1
          fprintf(fp, "%d %d -%.4f -%d %d %d -%d %d imp:n=1\n",
            cell_start + key, cell_start + key,
            array[key - 1].rho,
            2 + surface_start + z_cnt,
            coeff2*(k + 1), coeff1*k,
            i + 1, i);
          array[key - 1].cell_id = cell_start + key;
          key++;
        }
        else if (j == 2 + surface_start + z_cnt && (kind != 1))
        {
          fprintf(fp, "%d %d -%.4f -%d %d %d %d -%d %d imp:n=1\n",
            cell_start + key, cell_start + key,
            array[key - 1].rho,
            j, outer_cylinder_backup,
            coeff2*(k + 1), coeff1*k,
            i + 1, i);
          array[key - 1].cell_id = cell_start + key;
          key++;
        }
        else
        {
          // cell_number material_number density R-DIR XY-DIR Z-DIR imp:n=1
          fprintf(fp, "%d %d -%.4f -%d %d %d %d -%d %d imp:n=1\n",
            cell_start + key, cell_start + key,
            array[key - 1].rho,
            j, j - 1,
            coeff2*(k + 1), coeff1*k,
            i + 1, i);
          array[key - 1].cell_id = cell_start + key;
          key++;
        }
      }
    }
  }
  // close the file
  fclose(fp);
}

// write the corner region surface card and the cell card
// surface card and cell card generation
int McnpCornerSurfaceCellInp(
  double *face_z, int z_cnt,
  int surface_start,
  int cell_start,
  double R_OUT,
  double x_offset,    // the offset on X-axis of the cylinder circle location
  double y_offset,    // the offset on Y-axis of the cylinder circle location
  double up_tgx, 
  int CornerXcnt,
  PPOINT array)
{
  int i, j, k;
  int surface_number = 1;
  int cell_number = 1;
  double dx;
  FILE *fp_surface = NULL;
  FILE *fp_cell = NULL;
  fp_surface = fopen("CornerRegionSurface.txt", "a");
  fp_cell = fopen("CornerRegionCell.txt", "a");
  if (fp_surface == NULL || fp_cell == NULL)
  {
    Message("[*]McnpSurfaceCardPartInp : can not open the file\n");
    return -1;
  }
  dx = up_tgx / (double)CornerXcnt;
  // the up-right corner region
  for (k = 0; k < z_cnt; k++)
  {
    for (j = 0; j < CornerXcnt; j++)
    {
      for (i = j; i < CornerXcnt; i++)
      {
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + dx*i);       // x(+)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + dx*(i + 1)); // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n", surface_start + surface_number,
          y_offset + CIRCLE_FUN_Y(R_OUT, dx*j));                  // y(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n", surface_start + surface_number,
          y_offset + CIRCLE_FUN_Y(R_OUT, dx*(j + 1)));            // y(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        if (i == j)
        {
          fprintf(fp_surface, "%d c/z %.4f %.4f %.12f\n",
            surface_start + surface_number, x_offset, y_offset, R_OUT); // r
          surface_number++;
          // write the cell card
          // description:                x   x  y   y  z   z  r
          fprintf(fp_cell, "%d %d -%.4f %d -%d -%d %d %d -%d %d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number - 1].rho, // density
            surface_start + surface_number - 7,
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
        else
        {
          // write the cell card
          fprintf(fp_cell, "%d %d -%.4f %d -%d -%d %d %d -%d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number - 1].rho, // density
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
      }
    }
  }
  //the up-left corner region
  for (k = 0; k < z_cnt; k++)
  {
    for (j = 0; j < CornerXcnt; j++)
    {
      for (i = j; i < CornerXcnt; i++)
      {
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + (-1.0)*dx*i);   // x(-)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + (-1.0)*dx*(i + 1)); // x(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + CIRCLE_FUN_Y(R_OUT, dx*j)); // y(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + CIRCLE_FUN_Y(R_OUT, dx*(j + 1))); // y(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        if (i == j)
        {
          fprintf(fp_surface, "%d c/z %.4f %.4f %.12f\n",
            surface_start + surface_number, x_offset, y_offset, R_OUT); // r
          surface_number++;
          // write the cell card
          //                            x   x  y   y  z   z  r
          fprintf(fp_cell, "%d %d -%.4f -%d %d -%d %d %d -%d %d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number-1].rho, // density
            surface_start + surface_number - 7,
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
        else
        {
          // write the cell card
          fprintf(fp_cell, "%d %d -%.4f -%d %d -%d %d %d -%d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number - 1].rho, // density
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
      }
    }
  }
  //the bottom-left corner region
  for (k = 0; k < z_cnt; k++)
  {
    for (j = 0; j < CornerXcnt; j++)
    {
      for (i = j; i < CornerXcnt; i++)
      {
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + (-1.0)*dx*i);   // x(-)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + (-1.0)*dx*(i + 1)); // x(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + (-1.0)*CIRCLE_FUN_Y(R_OUT, dx*j)); // y(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + (-1.0)*CIRCLE_FUN_Y(R_OUT, dx*(j + 1))); // y(-)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        if (i == j)
        {
          fprintf(fp_surface, "%d c/z %.4f %.4f %.12f\n",
            surface_start + surface_number, x_offset, y_offset, R_OUT);
          surface_number++;
          // write the cell card
          //                            x   x  y   y  z   z  r
          fprintf(fp_cell, "%d %d -%.4f -%d %d %d -%d %d -%d %d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number-1].rho, // density
            surface_start + surface_number - 7,
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
        else
        {
          // write the cell card
          fprintf(fp_cell, "%d %d -%.4f -%d %d %d -%d %d -%d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number-1].rho, // density
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
      }
    }
  }
  //the bottom-right corner region
  for (k = 0; k < z_cnt; k++)
  {
    for (j = 0; j < CornerXcnt; j++)
    {
      for (i = j; i < CornerXcnt; i++)
      {
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + dx*i);   // x(+)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + dx*(i + 1)); // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + (-1.0)*CIRCLE_FUN_Y(R_OUT, dx*j)); // y(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + (-1.0)*CIRCLE_FUN_Y(R_OUT, dx*(j + 1))); // y(-)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        if (i == j)
        {
          fprintf(fp_surface, "%d c/z %.4f %.4f %.12f\n",
            surface_start + surface_number, x_offset, y_offset, R_OUT);
          surface_number++;
          // write the cell card
          //                            x   x  y   y  z   z  r
          fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d %d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number-1].rho, // density
            surface_start + surface_number - 7,
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
        else
        {
          // write the cell card
          fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
            cell_start + cell_number, cell_start + cell_number,
            array[cell_number-1].rho, // density
            surface_start + surface_number - 6,
            surface_start + surface_number - 5,
            surface_start + surface_number - 4,
            surface_start + surface_number - 3,
            surface_start + surface_number - 2,
            surface_start + surface_number - 1);
          array[cell_number - 1].cell_id = cell_start + cell_number;
          cell_number++;
        }
      }
    }
  }
  fclose(fp_surface);
  fclose(fp_cell);
  return surface_number;
}


// the four rectangular bar MCNP surface card and cell card generation
int McnpRectSurfaceCellInp(
  double *face_z, 
  int z_cnt,
  int a_cnt, int b_cnt,
  int surface_start,
  int cell_start,
  double x_offset,
  double y_offset,
  PPOINT array)
{
  int i, j, k;
  int surface_number = 1;
  int cell_number = 1;
  double da;
  double db;
  FILE *fp_surface = NULL;
  FILE *fp_cell = NULL;
  fp_surface = fopen("RectangularBarSurface.txt", "a");
  fp_cell = fopen("RectangularBarCell.txt", "a");
  if (fp_surface == NULL || fp_cell == NULL){
    Message("[*] McnpRectSurfaceCellInp : can not open the file...\n");
    return -1;
  }
  da = (OUTER_X_UP - BOARD_X_UP) / (double)a_cnt;
  db = (BOARD_X_UP - BOARD_X_DOWN) / (double)b_cnt;

  // first is the left rect region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < b_cnt; j++){
      for (i = 0; i < a_cnt; i++){
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_UP + da*i);          // x(+)
        surface_number++;
        if (i == a_cnt - 1){
          fprintf(fp_surface, "*%d px %.12f\n",
            surface_start + surface_number, x_offset + BOARD_X_UP + da*(i + 1));    // x(-)
        }
        else{
          fprintf(fp_surface, "%d px %.12f\n",
            surface_start + surface_number, x_offset + BOARD_X_UP + da*(i + 1));    // x(-)
        }
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_DOWN + db*j);        // y(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_DOWN + db*(j + 1));  // y(-)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);                             // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);                         // z(-)
        surface_number++;
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho,
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  // second is the top rect region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < a_cnt; j++){
      for (i = 0; i < b_cnt; i++){
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_DOWN + db*i);          // x(+)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_DOWN + db*(i + 1));    // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + da*j);            // y(+)
        surface_number++;
        if (j == a_cnt - 1){
          fprintf(fp_surface, "*%d py %.12f\n",
            surface_start + surface_number, y_offset + BOARD_Y_UP + da*(j + 1));      // y(-)
        }
        else{
          fprintf(fp_surface, "%d py %.12f\n",
            surface_start + surface_number, y_offset + BOARD_Y_UP + da*(j + 1));      // y(-)
        }
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);                               // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);                           // z(-)
        surface_number++;
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho,
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  // third is the left rect region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < b_cnt; j++){
      for (i = 0; i < a_cnt; i++){
        // write the surface card
        if (i == 0){
          fprintf(fp_surface, "*%d px %.12f\n",
            surface_start + surface_number, x_offset + OUTER_X_DOWN + da*i);        // x(+)
        }
        else{
          fprintf(fp_surface, "%d px %.12f\n",
            surface_start + surface_number, x_offset + OUTER_X_DOWN + da*i);        // x(+)
        }
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + da*(i + 1));  // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_DOWN + db*j);        // y(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_DOWN + db*(j + 1));  // y(-)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);                             // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);                         // z(-)
        surface_number++;
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho,
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  // fourth is the bottom rect region
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < a_cnt; j++){
      for (i = 0; i < b_cnt; i++){
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_DOWN + db*i);          // x(+)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_DOWN + db*(i + 1));    // x(-)
        surface_number++;
        if (j == 0){
          fprintf(fp_surface, "*%d py %.12f\n",
            surface_start + surface_number, y_offset + OUTER_Y_DOWN + da*j);          // y(+)
        }
        else{
          fprintf(fp_surface, "%d py %.12f\n",
            surface_start + surface_number, y_offset + OUTER_Y_DOWN + da*j);          // y(+)
        }
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + OUTER_Y_DOWN + da*(j + 1));    // y(-)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);                               // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);                           // z(-)
        surface_number++;
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho,
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  fclose(fp_surface);
  fclose(fp_cell);
  return surface_number;
}


// the four rectangular corner MCNP surface card and cell card generation
void McnpRectCornerSurfaceCellInp(
  double *face_z,
  int z_cnt,
  int cnt,
  int surface_start,
  int cell_start,
  double x_offset, double y_offset,
  PPOINT array)
{
  int i, j, k;
  double dx, dy;

  int surface_number = 1;
  int cell_number = 1;

  FILE *fp_surface = NULL;
  FILE *fp_cell = NULL;
  fp_surface = fopen("RectCornerSurface.txt", "a");
  fp_cell = fopen("RectCornerCell.txt", "a");
  
  if (fp_surface == NULL || fp_cell == NULL)
  {
    Message("[*] McnpRectCornerSurfaceCellInp : can not open the file...\n");
    return;
  }

  dx = dy = (OUTER_X_UP - BOARD_X_UP) / (double)cnt;
  
  // first is the up right corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_UP + dx*i); // x(+)
        surface_number++;
        if (i == cnt - 1){  // symmetry boundary
          fprintf(fp_surface, "*%d px %.12f\n",
            surface_start + surface_number, x_offset + BOARD_X_UP + dx*(i + 1)); // x(-)
        }else{
          fprintf(fp_surface, "%d px %.12f\n",
            surface_start + surface_number, x_offset + BOARD_X_UP + dx*(i + 1)); // x(-)
        }
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + dy*j);  // y(+)
        surface_number++;
        if (j == cnt - 1){ // symmetry boundary
          fprintf(fp_surface, "*%d py %.12f\n",
            surface_start + surface_number, y_offset + BOARD_Y_UP + dy*(j + 1)); // y(-)
        }else{
          fprintf(fp_surface, "%d py %.12f\n",
            surface_start + surface_number, y_offset + BOARD_Y_UP + dy*(j + 1)); // y(-)
        }
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        // write the cell card
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho, // density
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  // second is the up left corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        // write the surface card
        if (i == 0){
          fprintf(fp_surface, "*%d px %.12f\n",
            surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*i); // x(+)
        }
        else{
          fprintf(fp_surface, "%d px %.12f\n",
            surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*i); // x(+)
        }
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*(i + 1)); // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + dy*j);  // y(+)
        surface_number++;
        if (j == cnt - 1){
          fprintf(fp_surface, "*%d py %.12f\n",
            surface_start + surface_number, y_offset + BOARD_Y_UP + dy*(j + 1)); // y(-)
        }
        else{
          fprintf(fp_surface, "%d py %.12f\n",
            surface_start + surface_number, y_offset + BOARD_Y_UP + dy*(j + 1)); // y(-)
        }
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        // write the cell card
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho, // density
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  // third is the down left corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        // write the surface card
        if (i == 0){
          fprintf(fp_surface, "*%d px %.12f\n",
            surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*i); // x(+)
        }
        else{
          fprintf(fp_surface, "%d px %.12f\n",
            surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*i); // x(+)
        }
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*(i + 1)); // x(-)
        surface_number++;
        if (j == 0){
          fprintf(fp_surface, "*%d py %.12f\n",
            surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*j);  // y(+)
        }
        else{
          fprintf(fp_surface, "%d py %.12f\n",
            surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*j);  // y(+)
        }
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*(j + 1)); // y(-)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        // write the cell card
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho, // density
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  // fourth is the down right corner
  for (k = 0; k < z_cnt; k++){
    for (j = 0; j < cnt; j++){
      for (i = 0; i < cnt; i++){
        // write the surface card
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_UP + dx*i); // x(+)
        surface_number++;
        if (i == cnt - 1){
          fprintf(fp_surface, "*%d px %.12f\n",
            surface_start + surface_number, x_offset + BOARD_X_UP + dx*(i + 1)); // x(-)
        }
        else{
          fprintf(fp_surface, "%d px %.12f\n",
            surface_start + surface_number, x_offset + BOARD_X_UP + dx*(i + 1)); // x(-)
        }
        surface_number++;
        if (j == 0){
          fprintf(fp_surface, "*%d py %.12f\n",
            surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*j);  // y(+)
        }
        else{
          fprintf(fp_surface, "%d py %.12f\n",
            surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*j);  // y(+)
        }
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*(j + 1)); // y(-)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k]);      // z(+)
        surface_number++;
        fprintf(fp_surface, "%d pz %.12f\n",
          surface_start + surface_number, face_z[k + 1]);  // z(-)
        surface_number++;
        // write the cell card
        // write the cell card
        fprintf(fp_cell, "%d %d -%.4f %d -%d %d -%d %d -%d imp:n=1\n",
          cell_start + cell_number, cell_start + cell_number,
          array[cell_number - 1].rho, // density
          surface_start + surface_number - 6,
          surface_start + surface_number - 5,
          surface_start + surface_number - 4,
          surface_start + surface_number - 3,
          surface_start + surface_number - 2,
          surface_start + surface_number - 1);
        array[cell_number - 1].cell_id = cell_start + cell_number;
        cell_number++;
      }
    }
  }
  fclose(fp_surface);
  fclose(fp_cell);
}

// get the MCNP fuel database
// from 300K to 1800K, 25K per database
void McnpFuelDatabase(PPOINT array, int total)
{
  int i;
  // database number
  int data_base[62];
  int temp;
  double ratio;

  // make the data number
  for(i=0; i<62; i++)
  {
    data_base[i] = 10 + i;
  }
  Message("[*] ###### DEBUG temperature=%.4f\n", array[3].u.temperature);

  // begin to generate the data base number
  for(i=0; i<total; i++){
    if(array[i].u.temperature <= 300.0){
      ratio = 1.0;
      array[i].database = data_base[0];
    }else if(array[i].u.temperature >= 1800.0){
      ratio = 1.0;
      array[i].database = data_base[61];
    }else{
      temp = (int)(array[i].u.temperature - 300.0);
      temp = temp / 25;
      // use the formula to calculate the ratio
      ratio = (pow((325.0 + 25.0*temp), 0.5) - pow(array[i].u.temperature, 0.5)) /
              (pow((325.0 + 25.0*temp), 0.5) - pow((300.0 + 25.0*temp), 0.5));
      // choose the high temperature data base or low temperature data base
      if (0.5 < ratio && ratio <= 1.0){
        array[i].database = data_base[temp];   // low temperature data base
      }else{
        array[i].database = data_base[temp + 1];  // high temperature data base
      }
    }
  }
}

// generate the MCNP coolant database
// from 300K to 900K, 25K per database
void McnpWaterDatabase(PPOINT array, unsigned int total)
{
  int i;
  int database[26];
  double temp;

  for(i=0;i<26;i++){
    database[i] = 10+i;
  }
  for(i=0;i<total;i++){
    temp = array[i].u.temperature;
    if(temp<=300.0){
      array[i].database = database[0];
    }else if(temp >= 900.0){
      array[i].database = database[25];
    }else{
      temp = (temp - 300.0)/25.0;
      array[i].database = database[(int)temp];
    }
  }
  Message("[*] McnpWaterDatabase has done...\n");
}


// generate the MCNP material card
void McnpMaterialCardPartInp(const char *filename, 
  PPOINT array, 
  int total,
  int kind)
{
  unsigned int i;
  FILE *fp = fopen(filename, "a");
  if(fp == NULL)
  {
    Message("[*] McnpMaterialCardPartInp : can not open the file...\n");
    return ;
  }
  // determine which part material card will be written.
  switch (kind)
  {
    case WATER_TYPE:
    {
      for(i=0; i<total; i++)
      {
        fprintf(fp, "m%d   8016.%dc %f & \n", array[i].cell_id, array[i].database, 1.0);
        fprintf(fp, "      1001.%dc %f\n", array[i].database, 2.0);
      }
      break;
    }
    case FUEL_TYPE:
    {
      for(i=0; i<total; i++)
      {
        fprintf(fp, "m%d   8016.%dc 1.98 &\n", array[i].cell_id, array[i].database);
        fprintf(fp, "      92234.%dc 0.00055 & \n", array[i].database);
        fprintf(fp, "      92235.%dc 0.05 &\n", array[i].database);
        fprintf(fp, "      92238.%dc 0.949945\n", array[i].database);
      }
      break;
    }
    case CLAD_TYPE: // 900K
    {
      for(i=0; i<total; i++)
      {
        fprintf(fp, "m%d   24050.73c 4.3450e-05 \n", array[i].cell_id);
        fprintf(fp,"       24052.73c 8.3789e-04\n");
        fprintf(fp,"       24053.73c 9.5010e-05\n");
        fprintf(fp,"       24054.73c 2.3650e-05\n");
        fprintf(fp,"       26054.73c 1.1600e-04\n");
        fprintf(fp,"       26056.73c 1.8344e-03\n");
        fprintf(fp,"       26057.73c 4.4000e-05\n");
        fprintf(fp,"       26058.73c 5.6000e-06\n");
        fprintf(fp,"       50112.73c 1.4550e-04\n");
        fprintf(fp,"       50114.73c 9.9000e-05\n");
        fprintf(fp,"       50115.73c 5.1000e-05\n");
        fprintf(fp,"       50116.73c 2.1810e-03\n");
        fprintf(fp,"       50117.73c 1.1520e-03\n");
        fprintf(fp,"       50118.73c 3.6330e-03\n");
        fprintf(fp,"       50119.73c 1.2885e-03\n");
        fprintf(fp,"       50120.73c 4.8870e-03\n");
        fprintf(fp,"       50122.73c 6.9450e-04\n");
        fprintf(fp,"       50124.73c 8.6850e-04\n");
        fprintf(fp,"       40090.73c 0.5052390 \n");
        fprintf(fp,"       40091.73c 0.1101804 \n");
        fprintf(fp,"       40092.73c 0.1684130 \n");
        fprintf(fp,"       40094.73c 0.1706716 \n");
        fprintf(fp,"       40096.73c 0.0274960 \n");
      }
      break;
    }
    default:
      Message("Please recheck the kind number....\n");
  }
  // close the file
  fclose(fp);
  Message("[*] McnpMaterialCardPartInp has done...\n");
}


// This function is objective to generate the F7, SD7, KCODE etc. of the current part for MCNP data card input.
// PARAMETER : const char *filename => the filename of the F7 part card, which will be written with "a+".
//             PPOINT array => the structure array of the current points.
void McnpDataCardPartInp(const char *filename, int fuel_total)
{
  int i;
  int index = 1;
  int cnt = fuel_total;
  FILE *fp = fopen(filename, "a+");
  if (fp == NULL)
  {
    Message("[*] McnpF7CardPartInp : Open file error...\n");
    exit(0);
  }

  // judge the current fuel part's cell number whether can be divided by 5
  switch (cnt % 5)
  {
    case 0:  // the cell number can be divided by 5
    {
      // write the "F7" header
      fprintf(fp, "F7:N %d %d %d %d %d &\n     ", FUEL_START, FUEL_START + 1, FUEL_START + 2, FUEL_START + 3, FUEL_START + 4);
      // write the cell number of the fuel region, one line has five elements
      for (i = FUEL_START + 5; i <= FUEL_TOTAL - 5; i++)
      {
        fprintf(fp, "%d ", i);
        if (index % 5 == 0)
        {
          fputs("&\n", fp);
          fputs("     ", fp);
        }
        index++;
      }
      // write the last five elements of the fuel cell number
      for (i = FUEL_TOTAL - 4; i <= FUEL_TOTAL; i++)
      {
        fprintf(fp, "%d ", i);
      }
      fputs("\n", fp);

      // write the "SD7" header
      fputs("SD7  1 1 1 1 1 &\n", fp);
      for (i = 1; i < (cnt / 5) - 1; i++)
      {
        fputs("     1 1 1 1 1 &\n", fp);
      }
      fputs("     1 1 1 1 1", fp);
      fputs("\n", fp);
      break;
    }
    default: // the cell number can not be divided by 5
    {
      // write the "F7" header.
      fprintf(fp, "F7:N %d %d %d %d %d &\n    ", FUEL_START, FUEL_START + 1, FUEL_START + 2, FUEL_START + 3, FUEL_START + 4);
      // write the cell number of the fuel region, one line has five elements
      for (i = FUEL_START + 5; i <= FUEL_TOTAL - 5; i++)
      {
        fprintf(fp, "%d ", i);
        if (index % 5 == 0)
        {
          fputs("&\n", fp);
          fputs("    ", fp);
        }
        index++;
      }
      // write the some residual cell number of the fuel cell number
      for (i = FUEL_TOTAL - 4; i <= FUEL_TOTAL; i++)
      {
        fprintf(fp, "%d ", i);
      }
      fputs("\n", fp);

      // write the "SD7" header
      fputs("SD7  1 1 1 1 1 &\n", fp);
      for (i = 1; i < cnt / 5; i++)
      {
        fputs("     1 1 1 1 1 &\n", fp);
      }
      fputs("     ", fp);
      for (i = 0; i < cnt % 5; i++)
      {
        fputs("1 ", fp);
      }
      fputs("\n", fp);
    }
  }

  // write the KCODE card
  fputs("kcode 10000 1.0 10 60\n", fp);
  fputs("ksrc 0.00311 0.00317 10.00137\n", fp);
  // close the file
  fclose(fp);
}

// Generate MCNP every part card
// return value : a integer array, which contains the surface number and the cell number
// In the McnpPartGrid function, KIND=1, KIND=2 and KIND=3 are all be used.
// If use the MCNP input file generation, only KIND=1 and KIND=2 have be used,
// KIND=3 will be used in the other function.
void GenerateMcnpForCylinderRegion(double R_IN, double R_OUT, int layer_cnt,
  double begin_z, double end_z, int z_cnt,
  int N,     // the number of the sections, 360 / N.
  int kind,
  double up_tgx, double down_tgx,
  double up_tgy, double down_tgy,
  int surface_start,
  int cell_start,
  double x_offset,
  double y_offset,
  int type,   // judge the fuel type or the clad type
  int *RetVal)
{
  int outer_cylinder = 0;
  static int outer_cylinder_backup = 0;   // caution please: must be "static int" type!
  int x_plate;
  int y_plate;
  int surface_sum = 0;
  int total = 0;

  PPOINT point = NULL;
  double *face_r = NULL;
  double *face_z = NULL;

  // allocate the memory
  face_z = (double *)malloc(sizeof(double)*(z_cnt + 1));
  face_r = (double *)malloc(sizeof(double)*layer_cnt);
  memset(face_z, 0, sizeof(double)*(z_cnt + 1));
  memset(face_r, 0, sizeof(double)*layer_cnt);

  total = layer_cnt*z_cnt*N;
  point = (PPOINT)malloc(sizeof(POINT)*total);
  memset(point, 0, sizeof(POINT)*total);

  // density for cell card (DEBUG USE)
  int i;
  if (type == CLAD_TYPE){
    for(i=0;i<total;i++){
      point[i].rho = 6.55;
    }
  }

  // generate the grid
  McnpCylinderPartGrid(R_IN, R_OUT, layer_cnt,
    begin_z, end_z, z_cnt,
    N,
    kind,
    face_r,
    face_z,
    up_tgx,
    down_tgx,
    up_tgy,
    down_tgy,
    point,
    x_offset,
    y_offset);

  // get the global fuel cell
  // the main purpose is to record the fuel cell location and the MCNP cell volume.
  if(type == FUEL_TYPE)
  {
    memcpy(g_coarse_fuel+g_index, point, sizeof(POINT)*total);
    g_index += total;
  }
  
  // generate the surface card
  surface_sum = McnpCylinderSurfaceCardInp("SurfacePart.txt",
    face_r, layer_cnt,
    face_z, z_cnt,
    N,
    kind,
    up_tgx,
    down_tgx,
    up_tgy,
    down_tgy,
    surface_start,
    &outer_cylinder,
    &x_plate,
    &y_plate,
    x_offset,
    y_offset);

  // IDW interpolation here
  if(type == FUEL_TYPE){
    Message("Begin to IDW_TO_MCNP, fuel region...\n");
    IDW_TO_MCNP(FUEL_THREAD_ID, point, total);
    Message("IDW_TO_MCNP has done...\n");
  }
  else{
    ; // clad do not need the IDW interpolation;
  }

  // generate the cell card
  McnpCylinderCellCardInp("CellPart.txt",
    layer_cnt, z_cnt,
    N,
    kind,
    surface_start,
    cell_start,
    outer_cylinder,
    outer_cylinder_backup,
    x_plate,
    y_plate,
    surface_start + surface_sum,
    point);

  if(type == FUEL_TYPE)
  {
    // generate the MCNP database
    McnpFuelDatabase(point, total);
  }
  
  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", point, total, type);

  // free the allocated arrays
  free(point);
  free(face_r);
  free(face_z);
  point = NULL;
  face_r = NULL;
  face_z = NULL;

  // let "outer_cylinder_backup" to be "outer_cylinder" 's value
  outer_cylinder_backup = outer_cylinder;

  RetVal[0] = surface_sum;
  RetVal[1] = total;
}


// for KIND=3, corner region
void GenerateMcnpForCornerRegion(int surface_start, int cell_start,
  double R_OUT,
  double up_tgx,
  double begin_z, double end_z,
  int z_cnt,
  double x_offset,
  double y_offset,
  int CornerXcnt,
  int *RetVal)// return the surface number and the cell number
{
  int line_cell_cnt = CornerXcnt;
  int total = 0;
  double *face_z = NULL;
  PPOINT array = NULL;
  
  int surface_sum = 0;

  // get the whole number of corner region discretazation.
  while (line_cell_cnt)
  {
    total += line_cell_cnt--;
  }
  total = 4 * total * z_cnt;

  // allocate the memory
  face_z = (double *)malloc(sizeof(double)*(z_cnt + 1));
  memset(face_z, 0, sizeof(double)*(z_cnt + 1));
  array = (PPOINT)malloc(sizeof(POINT)*total);
  memset(array, 0, sizeof(POINT)*total);

  // density for cell card (DEBUG USE)
  /*
  int i;
  for (i = 0; i < total; i++)
  {
    array[i].rho = 1.0;
  }
  */

  // meshing
  McnpCornerPartGrid(R_OUT, 
    up_tgx, 
    face_z, 
    begin_z, end_z, z_cnt, 
    x_offset, y_offset, 
    CornerXcnt, 
    array);

  // IDW interpolation here
  IDW_TO_MCNP(WATER_THREAD_ID, array, total);
    
  // write in the MCNP input file.
  surface_sum = McnpCornerSurfaceCellInp(face_z, z_cnt, 
    surface_start, cell_start, 
    R_OUT, 
    x_offset, y_offset, 
    up_tgx, 
    CornerXcnt, 
    array);

  // generate the mcnp coolant database
  McnpWaterDatabase(array, total);

  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", array, total, WATER_TYPE);

  // free the allocated memory
  free(face_z);
  free(array);

  RetVal[0] = surface_sum;
  RetVal[1] = total;
}

// generate MCNP input for rectangular bar region
void GenerateMcnpForRectRegion(int surface_start, int cell_start,
  double begin_z, double end_z,
  int z_cnt,
  int a_cnt, int b_cnt,
  double x_offset,
  double y_offset,
  int *RetVal)
{
  double *face_z = NULL;
  PPOINT array = NULL;

  int surface_sum = 0;
  int total = 4 * a_cnt*b_cnt*z_cnt;

  // allocate the memory
  face_z = (double *)malloc(sizeof(double)*(z_cnt + 1));
  memset(face_z, 0, sizeof(double)*(z_cnt + 1));
  array = (PPOINT)malloc(sizeof(POINT)*total);
  memset(array, 0, sizeof(POINT)*total);

  // density for cell card (DEBUG USE)
  /*
  int i;
  for (i = 0; i < total; i++)
  {
    array[i].rho = 1.0;
  }
  */

  // meshing
  OuterMainRectPartGrid(a_cnt, b_cnt, face_z, begin_z, end_z, z_cnt,
    x_offset, y_offset,
    array);

  // IDW here
  IDW_TO_MCNP(WATER_THREAD_ID, array, total);

  // write in the MCNP input file
  surface_sum = McnpRectSurfaceCellInp(face_z, z_cnt, 
    a_cnt, b_cnt, 
    surface_start, cell_start,
    x_offset, y_offset,
    array);

  // generate the water database
  McnpWaterDatabase(array, total);

  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", array, total, WATER_TYPE);

  // free the allocated memory
  free(face_z);
  free(array);

  RetVal[0] = surface_sum;
  RetVal[1] = total;
}


// generate MCNP input for rectangular corner region
void GenerateMcnpForRectCornerRegion(int surface_start, int cell_start,
  double begin_z, double end_z,
  int z_cnt,
  int cnt,
  double x_offset,
  double y_offset)
{
  double *face_z = NULL;
  PPOINT array = NULL;
  int total = 4 * cnt*cnt*z_cnt;

  // allocate the memory
  face_z = (double *)malloc(sizeof(double)*(z_cnt + 1));
  memset(face_z, 0, sizeof(double)*(z_cnt + 1));
  array = (PPOINT)malloc(sizeof(POINT)*total);
  memset(array, 0, sizeof(POINT)*total);

  // density for cell card (DEBUG USE)
  /*
  int i;
  for (i = 0; i < total; i++)
  {
    array[i].rho = 1.0;
  }
  */

  // meshing
  OuterCornerRectPartGrid(cnt, face_z, begin_z, end_z, z_cnt,
    x_offset, y_offset,
    array);

  // IDW here
  IDW_TO_MCNP(WATER_THREAD_ID, array, total);

  // write in the MCNP input file
  McnpRectCornerSurfaceCellInp(face_z, z_cnt, cnt, surface_start, cell_start,
    x_offset, y_offset,
    array);

  // generate the MCNP water database
  McnpWaterDatabase(array, total);
  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", array, total, WATER_TYPE);

  // free the allocated memory
  free(face_z);
  free(array);
}


void Initial_MCNP_input()
{
  double begin_z = 0.0;
  double end_z = 20.0;
  int z_cnt = 20;
  int before_surface_cnt = 0;
  int before_cell_cnt = 0;
  int layer_cnt;
  int N;
  int kind;
  double R_IN;
  double R_OUT;
  double up_tgx = 0.6;
  double down_tgx = -0.6;
  double up_tgy = 0.6;
  double down_tgy = -0.6;
  double x_offset = 0.0;
  double y_offset = 0.0;
  int RetVal[2] = { 0 };

  // make the g_index return to the zero.
  g_index = 0;
  int total = 0;

  /************************************************************************/
  /*                          TEST CODE                                   */
  /************************************************************************/
  
  /* FUEL REGION */
  // 1. cylinder with the circle central axis
  R_IN = 0.0;
  R_OUT = 0.1;
  layer_cnt = 2;
  N = 4;
  kind = 1;
  GenerateMcnpForCylinderRegion(R_IN, R_OUT, layer_cnt, begin_z, end_z, z_cnt, N, kind,
    up_tgx, down_tgx, up_tgy, down_tgy,
    before_surface_cnt,
    before_cell_cnt,
    x_offset,
    y_offset,
    FUEL_TYPE,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;    // caution : here "RetVal[0] - 1" not the "RetVal[0]"
  before_cell_cnt += RetVal[1];
  Message("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  
  // 2. cylinder annulus
  R_IN = 0.1;
  R_OUT = 0.3;
  layer_cnt = 8;
  N = 16;
  kind = 2;
  GenerateMcnpForCylinderRegion(R_IN, R_OUT, layer_cnt, begin_z, end_z, z_cnt, N, kind,
    up_tgx, down_tgx, up_tgy, down_tgy,
    before_surface_cnt,
    before_cell_cnt, 
    x_offset,
    y_offset,
    FUEL_TYPE,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;    // caution : here "RetVal[0] - 1" not the "RetVal[0]"
  before_cell_cnt += RetVal[1];
  Message("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);

  // 3. cylinder annulus
  R_IN = 0.3;
  R_OUT = 0.5;
  layer_cnt = 8;
  N = 32;
  kind = 2;
  GenerateMcnpForCylinderRegion(R_IN, R_OUT, layer_cnt, begin_z, end_z, z_cnt, N, kind,
    up_tgx, down_tgx, up_tgy, down_tgy,
    before_surface_cnt,
    before_cell_cnt,
    x_offset,
    y_offset,
    FUEL_TYPE,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;    // caution : here "RetVal[0] - 1" not the "RetVal[0]"
  before_cell_cnt += RetVal[1];
  Message("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);

  /* CLADDING REGION */
  // 4. cylinder annulus
  R_IN = 0.5;
  R_OUT = 0.6;
  layer_cnt = 2;
  N = 4;
  kind = 3;    // legacy parameter
  GenerateMcnpForCylinderRegion(R_IN, R_OUT, layer_cnt, begin_z, end_z, z_cnt, N, kind,
    up_tgx, down_tgx, up_tgy, down_tgy,
    before_surface_cnt,
    before_cell_cnt,
    x_offset,
    y_offset,
    CLAD_TYPE,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;
  before_cell_cnt += RetVal[1];
  Message("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  Message("Cylinder region has done...\n");

  /* COOLANT REGION */
  // 5. four corner region
  int CornerXcnt = 10;
  GenerateMcnpForCornerRegion(before_surface_cnt, before_cell_cnt, 
    R_OUT, 
    up_tgx, 
    begin_z, end_z, z_cnt, 
    x_offset, y_offset, 
    CornerXcnt,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;
  before_cell_cnt += RetVal[1];
  Message("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  Message("Four corner region has done...\n");


  // 6. Rectangular bar region
  int a_cnt = 5;
  int b_cnt = 20;
  GenerateMcnpForRectRegion(before_surface_cnt, before_cell_cnt,
    begin_z, end_z, z_cnt,
    a_cnt, b_cnt,
    x_offset, y_offset,
    RetVal);
  before_surface_cnt += RetVal[0] - 1; 
  before_cell_cnt += RetVal[1];
  Message("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  Message("#################### %d, %d\n", RetVal[0], RetVal[1]);
  Message("Four Rectangular bar has done...\n");

  // 7. Rectangular corner region
  int cnt = 2;
  GenerateMcnpForRectCornerRegion(before_surface_cnt, before_cell_cnt, 
    begin_z, end_z, z_cnt, 
    cnt, 
    x_offset, y_offset);
  Message("Four Rectangular corner has done...\n");

  // generate the F7, SD7, KCODE, KSRC data card
  McnpDataCardPartInp("DataCardPart.txt", FUEL_TOTAL);
  Message("KCODE data card has done...\n");
}

// combine files
void CombineFile(char *CylinderCellFile, char *CornerCellFile,
  char *RectBarCellFile, char *RectCornerCellFile,
  char *CylinderSurfaceFile, char *CornerSurfaceFile,
  char *RectBarSurfaceFile, char *RectCornerSurfaceFile,
  char *MaterialFile, char *DataCard)
{
  FILE *f_cylindercell      = NULL;
  FILE *f_cornercell        = NULL;
  FILE *f_rectbarcell       = NULL;
  FILE *f_rectcornercell    = NULL;
  FILE *f_cylindersurface   = NULL;
  FILE *f_cornersurface     = NULL;
  FILE *f_rectbarsurface    = NULL;
  FILE *f_rectcornersurface = NULL;
  FILE *f_material          = NULL;
  FILE *f_data              = NULL;
  FILE *fp_combine          = NULL;
  char buffer[80];

  f_cylindercell      = fopen(CylinderCellFile,      "r");
  f_cornercell        = fopen(CornerCellFile,        "r");
  f_rectbarcell       = fopen(RectBarCellFile,       "r");
  f_rectcornercell    = fopen(RectCornerCellFile,    "r");
  f_cylindersurface   = fopen(CylinderSurfaceFile,   "r");
  f_cornersurface     = fopen(CornerSurfaceFile,     "r");
  f_rectbarsurface    = fopen(RectBarSurfaceFile,    "r");
  f_rectcornersurface = fopen(RectCornerSurfaceFile, "r");
  f_material          = fopen(MaterialFile,          "r");
  f_data              = fopen(DataCard,              "r");
  fp_combine          = fopen("mcin", "w");

  if (f_cylindercell    == NULL ||
    f_cornercell        == NULL ||
    f_rectbarcell       == NULL ||
    f_rectcornercell    == NULL ||
    f_cylindersurface   == NULL ||
    f_cornersurface     == NULL ||
    f_rectbarsurface    == NULL ||
    f_rectcornersurface == NULL ||
    f_material          == NULL ||
    f_data              == NULL ||
    fp_combine          == NULL)
  {
    Message("[*] CombineFile : can not open the file...\n");
    return;
  }
  // write the title
  fputs("c Fluent MCNP coupled calculation\n", fp_combine);
  fputs("c Xi'an Jiaotong University - NuTHeL\n", fp_combine);
  // write the cell card
  fputs("c cell card\n", fp_combine);
  while(1)
  {
    fgets(buffer, 80, f_cylindercell);
    if(feof(f_cylindercell))break;
    fputs(buffer, fp_combine);
  }
  while(1)
  {
    fgets(buffer, 80, f_cornercell);
    if(feof(f_cornercell)) break;
    fputs(buffer, fp_combine);
  }
  while(1)
  {
    fgets(buffer, 80, f_rectbarcell);
    if(feof(f_rectbarcell)) break;
    fputs(buffer, fp_combine);
  }
  while(1)
  {
    fgets(buffer, 80, f_rectcornercell);
    if(feof(f_rectcornercell)) break;
    fputs(buffer, fp_combine);
  }
  fputs("99000 0 -99000:99001:-99002:99003:-99004:99005 imp:n=0\n", fp_combine);
  fputs("\n", fp_combine);
  // write the surface card
  fputs("c surface card\n", fp_combine);
  while(1)
  {
    fgets(buffer, 80, f_cylindersurface);
    if(feof(f_cylindersurface)) break;
    fputs(buffer, fp_combine);
  }
  while(1)
  {
    fgets(buffer, 80, f_cornersurface);
    if(feof(f_cornersurface)) break;
    fputs(buffer, fp_combine);
  }
  while(1)
  {
    fgets(buffer, 80, f_rectbarsurface);
    if(feof(f_rectbarsurface)) break;
    fputs(buffer, fp_combine);
  }
  while(1)
  {
    fgets(buffer, 80, f_rectcornersurface);
    if(feof(f_rectcornersurface)) break;
    fputs(buffer, fp_combine);
  }
  fprintf(fp_combine, "*%d px %.4f\n", 99000, OUTER_X_DOWN);
  fprintf(fp_combine, "*%d px %.4f\n", 99001, OUTER_X_UP);
  fprintf(fp_combine, "*%d py %.4f\n", 99002, OUTER_Y_DOWN);
  fprintf(fp_combine, "*%d py %.4f\n", 99003, OUTER_Y_UP);
  fprintf(fp_combine, "%d pz %.4f\n", 99004, OUTER_Z_DOWN);
  fprintf(fp_combine, "%d pz %.4f\n", 99005, OUTER_Z_UP);
  fputs("\n", fp_combine);
  // write the material card
  fputs("c material card\n", fp_combine);
  while (1)
  {
    fgets(buffer, 80, f_material);
    if (feof(f_material)) break;
    fputs(buffer, fp_combine);
  }
  // write the data card
  while (1)
  {
    fgets(buffer, 80, f_data);
    if(feof(f_data)) break;
    fputs(buffer, fp_combine);
  }
  fclose(fp_combine);
  fclose(f_cylindercell);
  fclose(f_cylindersurface);
  fclose(f_cornercell);
  fclose(f_cornersurface);
  fclose(f_rectbarcell);
  fclose(f_rectbarsurface);
  fclose(f_rectcornercell);
  fclose(f_rectcornersurface);
  fclose(f_material);
  fclose(f_data);
}


// FUNCTION : GetMcnpOutputPower, (warning : single plate)
//            which will extract the power density of each cell after the MCNP calculation.
// PARAMETER: 1.const char *filename => MCNP output file name
//            2.int fuel_flag => the start first part fuel cell id.
// actually in this function, the fuel_flag is 1.
void GetMcnpOutputPower(const char *filename, int fuel_flag)
{
  int i = 0;
  double power, residual;
  char read_string[80] = { 0 };
  char target[80] = { 0 };
  FILE *fp = fopen(filename, "r");
  FILE *result = fopen("power_result", "a");
  if (fp == NULL)
  {
    Message("[*] GetMcnpOutputPower : open file %s error...\n", filename);
    exit(0);
  }

  // set the target string, for example: "_cell__400__...." ( underline '_' is the example of space )
  sprintf(target, " cell  %d  ", fuel_flag);

  // begin to read the content of MCNP output file
  while (1)
  {
    // get the content
    fgets(read_string, 80, fp);
    if (feof(fp))
    {
      Message("[*] GetMcnpOutputPower : Read Done!\n");
      break;
    }

    // if detected the target string, then we need to get the power density
    if (strstr(read_string, target))
    {
      fgets(read_string, 80, fp);                  // get the newline "\n"
      fscanf(fp, "%lf %lf", &power, &residual);    // get the power density and the residual
      fgets(read_string, 80, fp);                  // get the newline "\r\n"
      g_coarse_fuel[0].u.power = power;
      for (i = 1; i < FUEL_TOTAL; i++)
      {
        fgets(read_string, 80, fp);              // get the "_cell__xxx...."
        fgets(read_string, 80, fp);              // get the newline "\r\n"
        fscanf(fp, "%lf %lf", &power, &residual);// get the power density and the residual
        fgets(read_string, 80, fp);              // get the newline "\r\n"
        g_coarse_fuel[i].u.power = power;
      }
      // we get all the contents, and we break the loop.
      break;
    }
  }
  // close the file
  fclose(fp);

  for(i=0;i<FUEL_TOTAL; i++)
  {
    fprintf(result, "%lf, %lf\n", g_coarse_fuel[i].u.power, g_coarse_fuel[i].mcnp_volume);
  }
  fclose(result);
  return;
}

// get the record the keff
// return value is keff.
double GetKeff(const char *filename, const char *keff_recorder)
{
  FILE *fp = NULL;
  FILE *fp_keff = NULL;   // used to record the keff in file.
  const char *target = "final keff = ";
  char read_string[80] = {0};
  char result[8] = {0};   // the string of keff is 7 characters, so we should set here char result[8];
  double keff = 0.0;
  char *location = NULL;
  int i = 0;
  // open file
  fp = fopen(filename, "r");
  if(fp == NULL)
  {
    Message("[*] GetKeff : open file error...\n");
    exit(0);
  }
  fp_keff = fopen(keff_recorder, "a");
  if(fp_keff == NULL)
  {
    Message("[*] GetKeff : Open keff recorder error...\n");
    exit(0);
  }

  // begin to read the file
  while(1)
  {
    fgets(read_string, 80, fp);
    if(feof(fp))
    {
      Message("[*] GetKeff : Read Done!\n");
      break;
    }
    // if detected the target string, then we need to get the keff
    if((location=strstr(read_string, target)))
    {
      for(i=0; i<7; i++)
      {
        result[i] = *(location + 14 + i);
      }
      keff = atof(result);  // transfer the string to double precision float number.
      break;
    }
  }

  // write the file
  fprintf(fp_keff, "%lf\n", keff);
  
  // close the file
  fclose(fp);
  fclose(fp_keff);
  return keff;
}

// re-arrange the MCNP output power
void McnpPowerArrange()
{
  double power_sum = 0.0;
  int i, j;
  double temp = 0.0;

  FILE *fp = fopen("real_power", "w");
  
  for (i=0; i<FUEL_TOTAL; i++)
  {
    power_sum = power_sum + g_coarse_fuel[i].u.power;
  }
  // get the real power density of the each cell.
  // transfer the power(W) to power density(W/m3).
  for (i=0; i<FUEL_TOTAL; i++)
  {
    g_coarse_fuel[i].u.power = 1.0E6 * POWER_TOTAL * g_coarse_fuel[i].u.power 
                              / power_sum / g_coarse_fuel[i].mcnp_volume; // (W/m3)
    fprintf(fp, "%lf\n", g_coarse_fuel[i].u.power);
  }
  fclose(fp);
}


// thermal-physical properties
DEFINE_SPECIFIC_HEAT(water_capacity,T,Tref,h,yi)
{
  double capacity2=0.0;
  // capacity2=0.2385*T*T-242.87*T+66453;
 capacity2=4.621885*pow(10.0,-5.0)*pow(T,4.0)-0.09820142*pow(T,3.0)+78.28443*pow(T,2.0)-27740.38*T+3.690159*pow(10,6.0);
   *h = capacity2 * (T-Tref);
  return capacity2;
}

DEFINE_PROPERTY(Waterdensity,c,t)
{
  double density2=745.0;
  double temper = C_T(c,t);
  if(temper<=277.6)
     density2=1010.0;
  else if(temper<=600.0&&temper>277.6) 
     density2=(-4.118*pow(temper,4)*pow(10.0,-8.0)+6.74*pow(temper,3)*pow(10.0,-5)-4.306*temper*temper/100.0+11.82*temper-154.4);
  else 
     density2=657.4;
  return density2;
}

// energy source term
DEFINE_SOURCE(fuel_power, c, t, dS, eqn)
{
  double source;
  source=C_UDMI(c,t,0);
  dS[eqn] = 0.0; 
  return source; 
}

// allocate the memory for global array
DEFINE_INIT(my_init, d)
{
  Message("[*] Hello, this is DEFINE_INIT!\n");
  InitializeGlobal();
  Message("[*] The MCNP global parameter has been initialized!\n");
}

// release the allocated memory
DEFINE_EXECUTE_AT_EXIT(garbage_collection)
{
  DestroyGlobal();
  Message("[*] The allocated memory has been freed successfully!\n");
}


DEFINE_EXECUTE_AT_END(my_at_end)
{
  Domain *d;
  d = Get_Domain(1);
  double keff[0xFF] = {0};
  int i = 0;
  Message("[#] DEBUG MESSAGE : THE RESIDUAL IS LARGER THAN THE STYLIZED RESIDUAL!\n");
  if(strlen(DOMAIN_EQN_LABEL(d, 0)) > 0)
  {
    // if the continuty residual is smaller than stylized residual, 
    // the we should begin the MCNP calculation.
    if(((d->rtrms[0])/(d->rscale[0])) < residual)
    {
      Message("[#####] DEBUG MESSAGE: THE RESIDUAL IS LITTLE THAN THE STYLIZED RESIDUAL!\n");
      Message("[#####] NOW WE BEGIN TO RUN THE MCNP PROGRAM!\n");
      // generate the MCNP input file
      Initial_MCNP_input();
      CombineFile("CellPart.txt", "CornerRegionCell.txt",
                  "RectangularBarCell.txt", "RectCornerCell.txt",
                  "SurfacePart.txt", "CornerRegionSurface.txt",
                  "RectangularBarSurface.txt",
                  "RectCornerSurface.txt",
                  "MaterialPart.txt", "DataCardPart.txt");
      // begin to run the MCNP
      system("run.bat");
      // get the fuel power of each cell of MCNP output file.
      GetMcnpOutputPower("mcino", FUEL_START);
      // get keff
      keff[i] = GetKeff("mcino", "keff");
      i++;
      if(i>=1){
        if(fabs(keff[i] - keff[i-1]) <= ERR){
          Message("[*] DEFINE_EXECUTE_AT_END : keff is converged, now we stop the calculation...\n");
          exit(0);
        }
      }
      // make the normalization of the MCNP output power 
      McnpPowerArrange();
      // use the IDW method to get the cell's source term in Fluent.
      Message("Begin the IDW_TO_FLUENT\n");
      IDW_TO_FLUENT(FUEL_THREAD_ID, g_coarse_fuel, FUEL_TOTAL);
      Message("IDW_TO_FLUENT has done...\n");
      // delete the temp *.txt files
      system("clean.bat");
      // reset the residual and continue the calculation.
      residual = 0.2 * (d->rtrms[0]/d->rscale[0]);
    }
  }
}