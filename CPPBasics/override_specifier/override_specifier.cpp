#include <iostream>

class A
{
public:
	virtual const char* getName() { return "A"; }
};

class B final : public A // обратите внимание на модификатор final здесь
{
public:
	virtual const char* getName() override { return "B"; }
};

class C : public B // ошибка компил€ции: нельз€ наследовать класс final
{
public:
	virtual const char* getName() override { return "C"; }
};

int main()
{
	return 0;
}