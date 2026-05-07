#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NS 30

void copy_string(char *input, char *output) {

  strcpy(output, input);

}

int main() {
  char input_string[] = "This is the string to copy.";
  // 15 was too small to hold the input string
  // raising the number is sufficient but this is a good
  // place to use #define
  char *output_string = (char *) malloc(NS*sizeof(char));

  copy_string(input_string, output_string);
  printf("%s %s\n", output_string, input_string);

  free(output_string);

  return 0;
}
