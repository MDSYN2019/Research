/* vsscanf.h
 */
#ifndef VSSCANF_H
#define VSSCANF_H

char* Advance(
          char* bufp  /* in */);

void vsscanf(
         char*    buf     /* in  */,
         char*    format  /* in  */,
         va_list  argp    /* out */);

#endif
