#include <iostream>

class A
{
public:
	virtual const char* getName() { return "A"; }
};

class B : public A
{
public:
	// Заметили final в конце? Это означает, что метод переопределить уже нельзя
	virtual const char* getName() override final { return "B"; } // всё хорошо, переопределение A::getName()
};

class C : public B
{
public:
	virtual const char* getName() override { return "C"; } // ошибка компиляции: переопределение метода B::getName(), который является final
};

int main()
{
	return 0;
}