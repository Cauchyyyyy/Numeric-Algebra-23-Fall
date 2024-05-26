#include <iostream>
#include <vector>
#include<Windows.h>
#include"Function.h"
#include"Exercise.h"
#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<iostream>
using namespace std;

int main() {
	for (int i = 2; i<21; i++)
	{
		cout << i << "-order Hilbert matrix's infinite norm condition number is:";
		exercise_1(i);
	}
	cout << endl;
	for (int i = 5; i < 31; i++) {
		exercise_2(i);
		cout << endl;
	}
	return 0;
}