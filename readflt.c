#include <stdio.h>

void main(void)
{
  float a;
  while(fread(&a, 1, sizeof(float), stdin))
    printf("%g\n",a);

}
