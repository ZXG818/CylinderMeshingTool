#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ND_ND 3
#define PI 3.1415926
#define ERR 0.0001
#define BOARD_X_UP    9.0
#define BOARD_X_DOWN -9.0
#define BOARD_Y_UP    9.0
#define BOARD_Y_DOWN -9.0
#define OUTER_X_UP    11.0
#define OUTER_X_DOWN -11.0
#define OUTER_Y_UP    11.0
#define OUTER_Y_DOWN -11.0

// MACRO FUNCTION
#define AD2DEG(theta) ((theta)*PI / 180.0)
#define CIRCLE_FUN_Y(R,x) (sqrt((R)*(R) - (x)*(x)))

// A Cylinder Mesh Method
typedef struct _POINT_
{
  int cell_id;
  int database;
  double rho;
  union
  {
    double power;
    double temperature;
  }u;
  double xyz[ND_ND];
}POINT, *PPOINT;


// global parameter
PPOINT g_coarse_fuel = NULL;

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
        if (j == 0 && kind == 1) { r = face_r[0] / 2.0; }
        // the circle which not contains the central point or the annulus
        else { r = (face_r[i] + face_r[i - 1]) / 2.0; }

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
    printf("[*] McnpSurfaceCardInp : open file error...\n");
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
    printf("[*] McnpCellPartCardInp : open file error...\n");
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
        if (theta_sum >= 0.0 && theta_sum <= 90.0 - d_theta)      { coeff1 = 1; coeff2 = -1; }
        else if (fabs(theta_sum - 90.0) <= ERR)                   { coeff1 = 1; coeff2 = 1; }
        else if (theta_sum > 90.0 && theta_sum <= 270.0 - d_theta){ coeff1 = -1; coeff2 = 1; }
        else if (fabs(theta_sum - 270.0) <= ERR)                  { coeff1 = -1; coeff2 = -1; }
        else if (theta_sum > 270.0 && theta_sum <= 360)           { coeff1 = 1; coeff2 = -1; }
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
    printf("[*]McnpSurfaceCardPartInp : can not open the file\n");
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
          fprintf(fp_cell, "%d %d %.4f -%d %d -%d %d %d -%d %d imp:n=1\n",
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
          fprintf(fp_cell, "%d %d %.4f -%d %d %d -%d %d -%d %d imp:n=1\n",
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
          fprintf(fp_cell, "%d %d %.4f %d -%d %d -%d %d -%d %d imp:n=1\n",
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
    printf("[*] McnpRectSurfaceCellInp : can not open the file...\n");
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
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_UP + da*(i + 1));    // x(-)
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
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + da*(j + 1));      // y(-)
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
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + da*i);        // x(+)
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
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + OUTER_Y_DOWN + da*j);          // y(+)
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
    printf("[*] McnpRectCornerSurfaceCellInp : can not open the file...\n");
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
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_UP + dx*(i + 1)); // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + dy*j);  // y(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + dy*(j + 1)); // y(-)
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
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*i); // x(+)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*(i + 1)); // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + dy*j);  // y(+)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + BOARD_Y_UP + dy*(j + 1)); // y(-)
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
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*i); // x(+)
        surface_number++;
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + OUTER_X_DOWN + dx*(i + 1)); // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*j);  // y(+)
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
        fprintf(fp_surface, "%d px %.12f\n",
          surface_start + surface_number, x_offset + BOARD_X_UP + dx*(i + 1)); // x(-)
        surface_number++;
        fprintf(fp_surface, "%d py %.12f\n",
          surface_start + surface_number, y_offset + OUTER_Y_DOWN + dy*j);  // y(+)
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

// TODO : Not completed.
void McnpMaterialCardPartInp(const char *filename, PPOINT array, int total)
{
  int i;
  FILE *fp = fopen(filename, "a");
  if (fp == NULL)
  {
    printf("[*] McnpMaterialCardPartInp : open file error...\n");
    return;
  }

  for (i = 0; i < total; i++)
  {
    fprintf(fp, "m%d 92235 %.4f\n", array[i].cell_id, 1.0);
  }
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
  for (i = 0; i < total; i++)
  {
    point[i].rho = 1.0;
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
  
  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", point, total);

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
  int i;
  for (i = 0; i < total; i++)
  {
    array[i].rho = 1.0;
  }

  // meshing
  McnpCornerPartGrid(R_OUT, 
    up_tgx, 
    face_z, 
    begin_z, end_z, z_cnt, 
    x_offset, y_offset, 
    CornerXcnt, 
    array);
  // write in the MCNP input file.
  surface_sum = McnpCornerSurfaceCellInp(face_z, z_cnt, 
    surface_start, cell_start, 
    R_OUT, 
    x_offset, y_offset, 
    up_tgx, 
    CornerXcnt, 
    array);

  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", array, total);

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
  int i;
  for (i = 0; i < total; i++)
  {
    array[i].rho = 1.0;
  }

  // meshing
  OuterMainRectPartGrid(a_cnt, b_cnt, face_z, begin_z, end_z, z_cnt,
    x_offset, y_offset,
    array);

  // write in the MCNP input file
  surface_sum = McnpRectSurfaceCellInp(face_z, z_cnt, 
    a_cnt, b_cnt, 
    surface_start, cell_start,
    x_offset, y_offset,
    array);

  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", array, total);

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
  int i;
  for (i = 0; i < total; i++)
  {
    array[i].rho = 1.0;
  }

  // meshing
  OuterCornerRectPartGrid(cnt, face_z, begin_z, end_z, z_cnt,
    x_offset, y_offset,
    array);

  // write in the MCNP input file
  McnpRectCornerSurfaceCellInp(face_z, z_cnt, cnt, surface_start, cell_start,
    x_offset, y_offset,
    array);

  // generate the material card
  McnpMaterialCardPartInp("MaterialPart.txt", array, total);

  // free the allocated memory
  free(face_z);
  free(array);
}


void Initial_MCNP_input()
{
  double begin_z = 0;
  double end_z = 10;
  int z_cnt = 2;
  int before_surface_cnt = 0;
  int before_cell_cnt = 0;
  int layer_cnt;
  int N;
  int kind;
  double R_IN;
  double R_OUT;
  double up_tgx = 9.0;
  double down_tgx = -9.0;
  double up_tgy = 9.0;
  double down_tgy = -9.0;
  double x_offset = 0.0;
  double y_offset = 0.0;
  int RetVal[2] = { 0 };
  FILE *fp = fopen("BoundaryFace.txt", "w");

  /************************************************************************/
  /*                          TEST CODE                                   */
  /************************************************************************/
  
  // 1. cylinder with the circle central axis
  R_IN = 0.0;
  R_OUT = 1.0;
  layer_cnt = 2;
  N = 4;
  kind = 1;
  GenerateMcnpForCylinderRegion(R_IN, R_OUT, layer_cnt, begin_z, end_z, z_cnt, N, kind,
    up_tgx, down_tgx, up_tgy, down_tgy,
    before_surface_cnt,
    before_cell_cnt,
    x_offset,
    y_offset,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;    // caution : here "RetVal[0] - 1" not the "RetVal[0]"
  before_cell_cnt += RetVal[1];
  printf("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  
  // 2. cylinder annulus
  R_IN = 1.0;
  R_OUT = 5.0;
  layer_cnt = 5;
  N = 8;
  kind = 2;
  GenerateMcnpForCylinderRegion(R_IN, R_OUT, layer_cnt, begin_z, end_z, z_cnt, N, kind,
    up_tgx, down_tgx, up_tgy, down_tgy,
    before_surface_cnt,
    before_cell_cnt, 
    x_offset,
    y_offset,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;    // caution : here "RetVal[0] - 1" not the "RetVal[0]"
  before_cell_cnt += RetVal[1];
  printf("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);

  // 3. cylinder annulus
  R_IN = 5.0;
  R_OUT = 9.0;
  layer_cnt = 8;
  N = 16;
  kind = 3;    // legacy parameter
  GenerateMcnpForCylinderRegion(R_IN, R_OUT, layer_cnt, begin_z, end_z, z_cnt, N, kind,
    up_tgx, down_tgx, up_tgy, down_tgy,
    before_surface_cnt,
    before_cell_cnt,
    x_offset,
    y_offset,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;
  before_cell_cnt += RetVal[1];
  printf("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  printf("Cylinder region has done...\n");

  // 4. four corner region
  int CornerXcnt = 5;
  GenerateMcnpForCornerRegion(before_surface_cnt, before_cell_cnt, 
    R_OUT, 
    up_tgx, 
    begin_z, end_z, z_cnt, 
    x_offset, y_offset, 
    CornerXcnt,
    RetVal);
  before_surface_cnt += RetVal[0] - 1;
  before_cell_cnt += RetVal[1];
  printf("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  printf("Four corner region has done...\n");


  // 5. Rectangular bar region
  int a_cnt = 2;
  int b_cnt = 5;
  GenerateMcnpForRectRegion(before_surface_cnt, before_cell_cnt,
    begin_z, end_z, z_cnt,
    a_cnt, b_cnt,
    x_offset, y_offset,
    RetVal);
  before_surface_cnt += RetVal[0] - 1; 
  before_cell_cnt += RetVal[1];
  printf("[DEBUG] %d, %d\n", before_surface_cnt, before_cell_cnt);
  printf("#################### %d, %d\n", RetVal[0], RetVal[1]);
  printf("Four Rectangular bar has done...\n");

  // 6. Rectangular corner region
  int cnt = 4;
  GenerateMcnpForRectCornerRegion(before_surface_cnt, before_cell_cnt, 
    begin_z, end_z, z_cnt, 
    cnt, 
    x_offset, y_offset);
  printf("Four Rectangular corner has done...\n");

  // before_surface_cnt += 2;
  // write the boundary surface
  fprintf(fp, "%d px %.4f\n", 99000, OUTER_X_DOWN);
  fprintf(fp, "%d px %.4f\n", 99001, OUTER_X_UP);
  fprintf(fp, "%d py %.4f\n", 99002, OUTER_Y_DOWN);
  fprintf(fp, "%d py %.4f\n", 99003, OUTER_Y_UP);
  // close the file
  fclose(fp);
}

int main(int argc, char *argv[])
{
  Initial_MCNP_input();
  printf("Done all work!\n");
}