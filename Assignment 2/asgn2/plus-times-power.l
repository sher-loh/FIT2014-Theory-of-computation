/*  flex  file for lexical analyser for PLUS-TIMES-POWER expressions,
    which involve simple arithmetic for nonnegative integers
    with addition, multiplication and powers.
    Graham Farr, Monash University
    Initially, copied from the Wikipedia page for  lex, then modified over time.
    Last updated:  11 September 2021
*/
    
/*** Definition section ***/

%{
/* C code to be copied verbatim */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void yyerror(char *);
#include "quat.h"
/*
#include "y.tab.h"
*/
%}


%%
    /*** Rules section ***/

            /* yytext is a string containing the matched text. */

    /* Nonnegative integers ...  */

    /*  For QUAT, this regexp needs to be replaced by one for nonnegative reals. */
(0|[1-9][0-9]*)   {
            printf("Token: NUMBER;  Lexeme: %s\n", yytext);
/*
            yylval.num = atof(yytext);
            return  NUMBER;
*/
        }


    /*  Next seven lines of code only relevant to QUAT and Problems 4-6.  */
(ReplaceThisTextByRegexpForNonnegativeIntegers)   {
            printf("Token: WHOLENUMBER;  Lexeme: %s\n", yytext);
/*
            yylval.iValue = atoi(yytext);
            return  WHOLENUMBER;
*/
        }


Power  {
            printf("Token: POWER;  Lexeme: %s\n", yytext);
/*
            yylval.str = strdup(yytext);
            return  POWER;
*/
         }



[+*(),]   {
            printf("Token and Lexeme: %s\n", yytext);
/*
            return  *yytext;
*/
        }

\n      {
            printf("Token and Lexeme: <newline>\n");
/*
            return  *yytext;
*/
        }


[ \t]   {    }     /*  skip whitespace  */


.       {
           printf("Invalid character.\n");
/*
           yyerror("invalid character");
*/
        }
     /*  will match any single character that does not match
         any of the above patterns
     */


%%
    /*** C Code section ***/

void yyerror(char *s) {
      fprintf(stderr, "line %d: yytext = %s.  Error msg: %s.\n", yylineno, yytext, s);
}


int yywrap(void) {
   return 1;
}


/*  comment out the function  main()  if using  lex  with  yacc   */
int main(void) {
    yylex();
    return 0;
}
