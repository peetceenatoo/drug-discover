#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#define OLD "Drugs.csv"
#define NEW "Drugs(clean).csv"

int main(){

  // Open the files
  FILE* old = fopen(OLD, "r");
  FILE* new = fopen(NEW, "w");
  if( ( old == NULL ) || ( new == NULL ) ) return -1;

  char temp;
  bool smile;

  // The flag smiles will be true until the first comma is reached, and then be false until the \n is reached


  smile = true;
  // Skip first row
  while( ( temp = getc(old) ) != '\n' );
  // Read all the smiles
  while( ( temp = getc(old) ) != EOF ){
    if( temp == ',' ) smile = false;
    if( temp == '\n' ) smile = true;
    if( smile ) fprintf(new, "%c", temp);
  }

  // Close the files
  fclose(old);
  fclose(new);


  return 0;
}