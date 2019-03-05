#include <stdio.h>
#include <stdlib.h>

#define BUFLEN 400

system_(str,len)
char str[];
long len;
{
    int i;
    char str2[BUFLEN];
    
    i = len - 1;
    for(;;){
        if ( i < 0 ){
            i = 0 ;
            break;
        }
        if ( str[i] != ' ' ) {
            i++;
            break;
        }
        i--;
    }

    if ( i >= BUFLEN ){
        fprintf(stderr,"Error in system\n");
        exit(1);
    }
    strncpy(str2,str,i);
    
    str2[i] = '\0';

    printf("system(\"%s\")\n",str2);
    system(str2);
}
