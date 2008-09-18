#include <stdio.h>
#include <stdlib.h>

#define DELTA 0.001

#define MAXLINE 1024
#define FORMAT2020 "%E  %E %E %E"

void compare(double a, double b, int line, char* label)
{
    double x = a > b ? a - b : b - a;
    if(x > DELTA)
    {
        printf("%s: %22.16E != %22.16E on line %d\n", label, a, b, line);
    }
}

int main(int argc, char** argv)
{
    int linecount;
    char *f1name, *f2name;
    FILE *f1, *f2;
    float x1, y1, z1, c1, x2, y2, z2, c2;
    char line[MAXLINE];
    
    
    /* Check arguments */
    if(argc < 3)
    {
        fprintf(stderr, "Usage: %s file1 file2\n", argv[0]);
        exit(1);
    }
    
    linecount = 0;
    f1name = argv[1];
    f2name = argv[2];
    
    if((f1 = fopen(f1name, "r")) == NULL)
    {
        fprintf(stderr, "Couldn't open file \"%s\".", f1name);
        exit(1);
    }
    
    if((f2 = fopen(f2name, "r")) == NULL)
    {
        fprintf(stderr, "Couldn't open file \"%s\".", f2name);
        exit(1);
    }
    
    while(1)
    {   
        ++linecount;
        
        if(fgets(line, MAXLINE, f1) == NULL) break;
        fscanf(f1, FORMAT2020, &x1, &y1, &z1, &c1);
        
        if(fgets(line, MAXLINE, f2) == NULL) break;
        fscanf(f2, FORMAT2020, &x2, &y2, &z2, &c2);
        
        compare(x1, x2, linecount, "x");
        compare(y1, y2, linecount, "y");
        compare(z1, z2, linecount, "z");
        compare(c1, c2, linecount, "c");
    }
        
    return 0;
}
