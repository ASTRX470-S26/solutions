#include <stdio.h>
#include <stdlib.h>

int* create_array(int size) {
  int* arr = (int*)malloc(size * sizeof(int));
  if (arr == NULL) {
    printf("Memory allocation failed\n");
    return NULL;
  }

  for (int i = 0; i < size; i++) {
    arr[i] = i;
  }

  return arr;
}

void print_array(int* arr, int size) {
  for (int i = 0; i < size; i++) {
    printf("%d ", arr[i]);
  }
  printf("\n");
}

void modify_array(int* arr, int size) {

  for (int i = 0; i < size; i++) {
    arr[i] = i*i;
  }
  free(arr);
}

int main() {
  
  int size = 10;
  int* my_array = create_array(size);
  
  if (my_array == NULL) {
    return 1;
  } 

  // print initial array
  print_array(my_array, size);

  // modify array
  modify_array(my_array, size);

  // print modified array
  print_array(my_array, size);
  
  return 0;
}
