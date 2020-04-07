#include <iostream>

class A
{
public:
	virtual const char* getName1(int x) { return "A"; }
	virtual const char* getName2(int x) { return "A"; }
	virtual const char* getName3(int x) { return "A"; }
};

class B : public A
{
public:
	virtual const char* getName1(short int x) override { return "B"; } // ошибка компил€ции, метод не €вл€етс€ переопределением
	virtual const char* getName2(int x) const override { return "B"; } // ошибка компил€ции, метод не €вл€етс€ переопределением
	virtual const char* getName3(int x) override { return "B"; } // всЄ хорошо, метод €вл€етс€ переопределением A::getName3(int)

};

int main()
{
	return 0;
}