#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void copy_string(char *input, char *output) {

  strcpy(output, input);

}

int main() {
  char input_string[] = "This is the string to copy.";
  char *output_string = (char *) malloc(15*sizeof(char));

  copy_string(input_string, output_string);
  printf("%s %s\n", output_string, input_string);

  free(output_string);

  return 0;
}
