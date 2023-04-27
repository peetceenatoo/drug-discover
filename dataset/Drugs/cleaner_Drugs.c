#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#define OLD "Drugs.csv"
#define NEW "Drugs(clean).csv"

int main(){
  FILE* old = fopen(OLD, "r");
  FILE* new = fopen(NEW, "w");
  if(old == NULL || new == NULL) return -1;
  char temp;
  bool smile;

  smile = true;
  while( (temp = getc(old)) != '\n');
  while( (temp = getc(old)) != EOF){
    if(temp == ',') smile=false;
    if(temp == '\n') smile=true;
    if(smile) fprintf(new, "%c", temp);
  }

  return 0;
}