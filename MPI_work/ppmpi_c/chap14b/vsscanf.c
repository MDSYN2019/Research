/* vsscanf.c
 *
 * Copies the contents of buf into the arguments in the
 *    variable length argument list.
 *
 * Only supports string, integer, float and double.  So
 *     format should only contain %s, %d, %f and %lf.
 */
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "vsscanf.h"

void vsscanf(
         char*    buf     /* in  */,
         char*    format  /* in  */,
         va_list  argp    /* out */) {

    char*    fmtp;
    char*    bufp;
    char*    sval_ptr;
    int*     ival_ptr;
    double*  dval_ptr;
    float*   fval_ptr;

#ifdef CDEBUG
    printf("In vsscanf\n");
#endif

    bufp = buf;

#ifdef CDEBUG
    printf("At start bufp = %s, format = %s\n", bufp, format);
#endif

    for (fmtp = format; *fmtp; fmtp++) {
        if (*fmtp == '%')
          switch (*++fmtp) {
            case 'd':
#ifdef CDEBUG
                printf("fmtp = %s, bufp = %s\n", fmtp, bufp);
#endif
                ival_ptr = va_arg(argp, int *);    /* advance argp */
                sscanf(bufp, "%d", ival_ptr);
                bufp = Advance(bufp);
                break;
            case 'f':
#ifdef CDEBUG
                printf("fmtp = %s, bufp = %s\n", fmtp, bufp);
#endif
                fval_ptr = va_arg(argp, float *);  /* advance argp */
                sscanf(bufp, "%f", fval_ptr);
                bufp = Advance(bufp);
                break;
            case 'l':
#ifdef CDEBUG
                printf("fmtp = %s, bufp = %s\n", fmtp, bufp);
#endif
                dval_ptr = va_arg(argp, double *);  /* advance argp */
                sscanf(bufp, "%lf", dval_ptr);
                bufp = Advance(bufp);
                break;
            case 's':
#ifdef CDEBUG
                printf("fmtp = %s, bufp = %s\n", fmtp, bufp);
#endif
                sval_ptr = va_arg(argp, char *);    /* advance argp */
                sscanf(bufp, "%s", sval_ptr);
                bufp = Advance(bufp);
                break;
            default:
#ifdef CDEBUG
                printf("fmtp = %s, bufp = %s\n", fmtp, bufp);
#endif
                break;
        }  /* switch */
    }  /* for */

}  /* vsscanf */


/********************************************************/
/* Advance buffer pointer to next valid character by
 *     1.  Skipping nonwhite characters.
 *     2.  Skipping white space.
 */

char* Advance(
          char* bufp /* in */) {

    char* new_bufp = bufp;

    /* Skip over nonwhite space */
    while ((*new_bufp != ' ')  && (*new_bufp != '\t') &&
           (*new_bufp != '\n') && (*new_bufp != '\0'))
        new_bufp++;

    /* Skip white space */
    while ((*new_bufp == ' ')  || (*new_bufp == '\t') ||
           (*new_bufp == '\n') || (*new_bufp == '\0'))
        new_bufp++;

    return new_bufp;
}  /* Advance */
