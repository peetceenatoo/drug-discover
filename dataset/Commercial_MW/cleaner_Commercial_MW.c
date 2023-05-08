#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>

#define OLD1 "Commercial_MWlower330.csv"
#define OLD2 "Commercial_MW330-500.csv"
#define OLD3 "Commercial_MWhigher500.csv"

#define NEW1 "Commercial_MWlower330(clean).csv"
#define NEW2_1 "Commercial_MW330-500(clean)1.csv"
#define NEW2_2 "Commercial_MW330-500(clean)2.csv"
#define NEW3 "Commercial_MWhigher500(clean).csv"

int main(){

  // Open the files
  FILE* old1 = fopen(OLD1, "r");
  FILE* old2 = fopen(OLD2, "r");
  FILE* old3 = fopen(OLD3, "r");
  FILE* new1 = fopen(NEW1, "w");
  FILE* new2_1 = fopen(NEW2_1, "w");
  FILE* new2_2 = fopen(NEW2_2, "w");
  FILE* new3 = fopen(NEW3, "w");
  if(old1 == NULL || old2 == NULL || old3 == NULL || new1 == NULL || new2_1 == NULL || new2_2 == NULL || new3 == NULL) return -1;

  char temp;
  bool smile;

  // The flag smiles will be true until the first comma is reached, and then be false until the \n is reached

  // File old1 -> new1
  smiles = true;
  // Skip first row
  while( ( temp = getc(old1) ) != '\n' );
  // Read all the smiles
  while( ( temp = getc(old1) ) != EOF ){
    if( temp == ',' ) smiles = false;
    else if( temp == '\n' ) smiles = true;
    if( smiles ) fprintf(new1, "%c", temp);
  }

  // Currently half of the rows in the dataset "Commercial_MW330-500.csv"
  // Need to calculate count runtime soon
  int count = 1407638;

  // File old2 -> new2_1 and new2_2
  smiles = true;
  // Skip first row
  while( ( temp = getc(old2) ) != '\n' );
  // Read all the smiles
  while( ( temp = getc(old2) ) != EOF ){
    if( temp == ',' ) smiles = false;
    else if( temp == '\n' ){
      smiles = true;
      count--;
      }
    // If half of the file is reached, print in new2_1, else in new2
    if( smiles && ( count>0 ) ) fprintf(new2_1, "%c", temp);
    else if( smiles && ( count<=0 ) ) fprintf(new2_2, "%c", temp);
  }

  // File old3 -> new3
  smiles = true;
  // Skip first row
  while( (temp = getc(old3)) != '\n');
  // Read all the smiles
  while( ( temp = getc(old3) ) != EOF ){
    if( temp == ',' ) smile = false;
    else if( temp == '\n' ) smile = true;
    if( smile ) fprintf(new3, "%c", temp);
  }

  // Close the files
  fclose(old1);
  fclose(old2);
  fclose(old3);
  fclose(new1);
  fclose(new2_1);
  fclose(new2_2);
  fclose(new3);

  return 0;
}