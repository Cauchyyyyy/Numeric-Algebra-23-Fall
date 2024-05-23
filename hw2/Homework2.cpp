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
		cout << i << "阶Hilbert矩阵的无穷范数条件数是";
		exercise_1(i);
	}
	cout << endl;
	for (int i = 5; i < 31; i++) {
		exercise_2(i);
		cout << endl;
	}
	return 0;
}