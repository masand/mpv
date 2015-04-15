#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>

#include "config.h"
#include "common/msg.h"

#include "video/img_format.h"
#include "video/mp_image.h"
#include "vf.h"


// the pointers to the lookup tables for the Y and UV plane
struct vf_priv_s {
    int *Y;
    int *UV;
};

static void fill_half_of_lookup_table(int width, int height,
                                      int **lookup_table, int left)
{
    // default: left side of the image -> center of that is 
    // width/4, height/2
    int xmid = width / 4;
    int ymid = height / 2;
    // start at upper left corner
    int startx = 0;
    int starty = 0;
    // stop at the middle of the lowest line
    int stopx = width / 2;
    int stopy = height;
    // not the left side -> start in the middle in x direction
    if (left != 1) {
        // right side -> center of that is 3*width/4, height/2
        xmid = 3 * width / 4;
        // start at the middle of the first line
        startx = width / 2;
        // stop at the lower right corner
        stopx = width;
    }
    
    for (int y = starty; y < stopy; y++) {
        for (int x = startx; x < stopx; x++) {
            /*
             *  NEXT STEP: Try with inverse function!
                y = 0.24x⁴ + 0.22x² + 1 for x > 0 (Radius!)
                x = 0.24y⁴ + 0.22y² + 1
                 substitute: t = y²
                x = 0.24t² + 0.22t + 1 = 0.24(t² + 0.22/0.24t + 1/0.24)
                  = 0.24(t² + 0.22/0.24t + (0.22/0.48)² - (0.22/0.48)² + 1/0.24)
                  = 0.24((t + 0.22/0.48)² - 0.22²/0.48² + 1/0.24)
                  = 0.24((t + 0.11/0.24)² - 0.0484/0.48² + 2*0.48/0.48²)
                  = 0.24((t + 0.11/0.24)² + 0.9116/0.48²)
                <=>
                x/0.24 - 0.9116/0.48² = (t + 0.11/0.24)²
                t = -0.11/0.24 +/- sqrt(x/0.24 - 0.9116/0.48²)
                =>
                y = +/- sqrt( -0.11/0.24 +/- sqrt(x/0.24 - 0.9116/0.48²))

                => in our case (only positive values possible):
                y = sqrt(-0.11/0.24 + sqrt(x/0.24 - 0.9116/0.48²))
                * 
                * Mathematica: -0.204124 Sqrt[-11. - 1. Sqrt[-2279. + 2400. #1]] &
             */
             /* Inverse 2: 
                y = 0.24x⁵ + 0.22x³ + x for x > 0 (Radius!)
                x = 0.24y⁵ + 0.22y³ + y
                 substitute: t = y²; y = sqrt(t)
                x = sqrt(t)*(0.24t² + 0.22t + 1) = sqrt(t)*(0.24(t² + 0.22/0.24t + 1/0.24))
                x/sqrt(t) = 0.24(t² + 0.22/0.24t + (0.22/0.48)² - (0.22/0.48)² + 1/0.24)
                          = 0.24((t + 0.22/0.48)² - 0.22²/0.48² + 1/0.24)
                          = 0.24((t + 0.11/0.24)² - 0.0484/0.48² + 2*0.48/0.48²)
                          = 0.24((t + 0.11/0.24)² + 0.9116/0.48²)
                <=>
                x/sqrt(t)/0.24 - 0.9116/0.48² = (t + 0.11/0.24)²
                t = -0.11/0.24 +/- (sqrt(x/0.24 - 0.9116/0.48²) / sqrt(y))
                =>
                y = +/- sqrt( -0.11/0.24 +/- sqrt(x/0.24 - 0.9116/0.48²))

                => in our case (only positive values possible):
                y = sqrt(-0.11/0.24 + sqrt(x/0.24 - 0.9116/0.48²))
                * 
                * Mathematica: -0.204124 Sqrt[-11. - 1. Sqrt[-2279. + 2400. #1]] &

              */
            // inspired by: http://jsfiddle.net/s175ozts/4/
            int myX = x-xmid;
            int myY = y-ymid;
            // according to the author on the website above, we need to
            // normalize the radius for proper calculation in the atan,
            // sin and cos functions. Also the distortion function needs
            // to work in [0,1]
            double maxR = sqrt(pow((stopx - startx) / 2.0, 2) +
                              pow(ymid, 2));
            // range of myR: [0, maxR]
            double myR = sqrt(pow(myX, 2) + pow(myY, 2));
            // range of myRN: [0,1]
            double myRN = myR / maxR;
            // use the distortion function
            // -> range of newR: [0, 1.46*maxR]
            double newR = myR * (0.24 * pow(myRN, 4) + 
                                 0.22 * pow(myRN, 2) + 1);
            // TODO: Do stop here if newR > maxR (worked - nearly...)
            /*if (newR > maxR) {
                // use -1 to mark as invalid
                (*lookup_table)[2 * y * width + 2 * x + 0] = -1;
                (*lookup_table)[2 * y * width + 2 * x + 1] = -1;
            }
            else {*/
            printf("%.10f -> %.10f\n", newR/maxR, myR/maxR);
            // then calculate back to x/y coordinates
            double alpha = atan2(myY, myX);
            double newXf = fabs(xmid + cos(alpha) * newR);
            double newYf = fabs(ymid + sin(alpha) * newR);
            // calculate the new radius to doublecheck (otherwise there
            // were some parts of the picture projected to the outside
            // the wished "barrel")
            double gnRadius = sqrt(pow(newXf - xmid, 2) + 
                                  pow(newYf - ymid, 2));
            int newX = (int) (newXf + 0.5);
            int newY = (int) (newYf + 0.5);
            // the above function gives us a double-concarve projection,
            // it is the inverse function of what we want but we use the
            // inverse to find out what would go outside the original
            // border. Those pixels need then to be set to black; the
            // real barrel projection is then created by swapping x and
            // newX when accessing the array values
            if(floor(newR) == floor(gnRadius) && newX >= startx &&
               newX < stopx && newY >= starty && newY < stopy)
            {
                (*lookup_table)[2 * y * width + 2 * x + 0] = newX;
                (*lookup_table)[2 * y * width + 2 * x + 1] = newY;
            } else {
                // use -1 to mark as invalid
                //(*lookup_table)[2 * y * width + 2 * x + 0] = -1;
                //(*lookup_table)[2 * y * width + 2 * x + 1] = -1;
            }
        }
        //}
    }
    ///// IDEA: Calculate x⁴ for x=0 to maxR, then build the inverse
    // values and check whether 0 to maxR returns!
    
    /*//for (int x = startx; x < stopx; x++) {
    for (int x = 0; x < 75; x++) {        
        //for (int y = starty; y < stopy; y++) {
        for (int y = 0; y < 75; y++) {
            int xlookup = (*lookup_table)[2 * y * width + 2 * x + 0];
            int ylookup = (*lookup_table)[2 * y * width + 2 * x + 1];
            if (xlookup == -1 || ylookup == -1)
                continue;
            printf("xlookup = %d, ylookup = %d\n", xlookup, ylookup);
            int myX = xlookup-xmid;
            int myY = ylookup-ymid;
            printf("myX = %d, myY = %d\n", myX, myY);
            double maxR = 1.46*(sqrt(pow((stopx - startx) / 2.0, 2) +
                               pow(ymid / 2.0, 2)));
            double myR = sqrt(pow(myX, 2) + pow(myY, 2));
            double myRN = myR / maxR;
            printf("myR = %f, myRN = %f\n", myR, myRN);
            // use the distortion function
            //sqrt(-0.11/0.24 + sqrt(x/0.24 - 0.9116/0.48²))
            double newR = myR * sqrt(-0.11/0.24 + 
                                    sqrt(myRN/0.24 - 0.9116/pow(0.48,2)));
            printf("newR = %f\n", newR);
            // then calculate back to x/y coordinates
            double alpha = atan2(myY, myX);
            double newXf = fabs(xmid + cos(alpha) * newR);
            double newYf = fabs(ymid + sin(alpha) * newR);
            // calculate the new radius to doublecheck (otherwise there
            // were some parts of the picture projected to the outside
            // the wished "barrel")
            double gnRadius = sqrt(pow(newXf - xmid, 2) + 
                                  pow(newYf - ymid, 2));
            int newX = (int) (newXf + 0.5);
            int newY = (int) (newYf + 0.5);
            if (newX != x) printf("newX = %d, x = %d\n", newX, x);
            if (newY != y) printf("newY = %d, y = %d\n", newY, y);
        }
    }*/
}

static void tmp_print_lookup_table(int width, int height, int **lookup_table)
{
    printf("P2\n%d %d\n255\n", width, height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int xlookup = (*lookup_table)[2 * y * width + 2 * x + 0];
            int ylookup = (*lookup_table)[2 * y * width + 2 * x + 1];
            if (xlookup == -1 || ylookup == -1) printf("0 ");
            else printf("255 ");
        }
        printf("\n");
    }
}

static int init_lookup_table(int width, int height, int **lookup_table)
{
    // we store x and y value -> "*2"
    *lookup_table = malloc(width * height * 2 * sizeof(int));
    if (!lookup_table)
        return 0;
    fill_half_of_lookup_table(width, height, lookup_table, 1);
    fill_half_of_lookup_table(width, height, lookup_table, 0);
    //tmp_print_lookup_table(width, height, lookup_table);
    return 1;
}

static int config(struct vf_instance *vf, int width, int height,
          int d_width, int d_height,
          unsigned int flags, unsigned int fmt)
{
    // In query_format, we accept only the IMGFMT_420P ->
    // we know we are working with YUV -> init both tables
    struct vf_priv_s *tables = vf->priv;
    if (!init_lookup_table(width, height, &(tables->Y)))
        return 0;
    if (!init_lookup_table(width / 2, height / 2, &(tables->UV)))
        return 0;
    return vf_next_config(vf, width, height, d_width, d_height, flags, fmt);
}

static void uninit(struct vf_instance *vf)
{
    struct vf_priv_s *tables = vf->priv;
    free(tables->Y);
    free(tables->UV);
}

static struct mp_image *filter(struct vf_instance *vf, struct mp_image *mpi)
{
    mp_image_t *dmpi = vf_alloc_out_image(vf);
    if (!dmpi)
        return NULL;
    mp_image_copy_attributes(dmpi, mpi);
    
    for (int p = 0; p < mpi->num_planes; p++) {
        struct vf_priv_s *tables = vf->priv;
        int* lookup_table = p ? tables->UV : tables->Y;
        int curPlaneWidth = mpi->plane_w[p];
        int curPlaneHeight = mpi->plane_h[p];
        for (int y = 0; y < curPlaneHeight; y++) {
            // in lookup_table, we stored the double concarve projection
            // but now we need the inverse, convex projection -> use
            // newX,newY in src image and copy that pixel to x,y in the
            // destination
            uint8_t *p_dst = dmpi->planes[p] + dmpi->stride[p] * y;         
            for (int x = 0; x < curPlaneWidth; x++) {
                int newX = *lookup_table++;
                int newY = *lookup_table++;
                uint8_t *p_src = mpi->planes[p] + mpi->stride[p] * newY;
                if (newX == -1) {
                    // set to "black" for Y and UV
                    p_dst[x] = (p==0) ? 0 : 128;
                } else {
                    p_dst[x] = p_src[newX];
                }
            }
        }
    }

    talloc_free(mpi);
    return dmpi;
}

static int query_format(struct vf_instance *vf, unsigned int fmt)
{
    if (fmt != IMGFMT_420P)
        return 0;
    return vf_next_query_format(vf, fmt);
}

static int vf_open(vf_instance_t *vf){
    vf->config=config;
    vf->filter=filter;
    vf->query_format=query_format;
    vf->uninit=uninit;
    return 1;
}

const vf_info_t vf_info_riftdk2 = {
    .description = "oculus rift dk2",
    .name = "riftdk2",
    .open = vf_open,
    .priv_size = sizeof(struct vf_priv_s),
};
