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

  smile = true;
  while( (temp = getc(old1)) != '\n');
  while( (temp = getc(old1)) != EOF){
    if(temp == ',') smile=false;
    if(temp == '\n') smile=true;
    if(smile) fprintf(new1, "%c", temp);
  }

  smile = true;
  while( (temp = getc(old2)) != '\n');
  int count = 1407638; //half of the file rows
  while( (temp = getc(old2)) != EOF){
    if(temp == ',') smile=false;
    if(temp == '\n') {
      smile=true;
      count--;
      }
    if(smile && count>0) fprintf(new2_1, "%c", temp);
    if(smile && count<=0) fprintf(new2_2, "%c", temp);
  }

  smile = true;
  while( (temp = getc(old3)) != '\n');
  while( (temp = getc(old3)) != EOF){
    if(temp == ',') smile=false;
    if(temp == '\n') smile=true;
    if(smile) fprintf(new3, "%c", temp);
  }

  fclose(old1);
  fclose(old2);
  fclose(old3);

  fclose(new1);
  fclose(new2_1);
  fclose(new2_2);
  fclose(new3);

  return 0;
}