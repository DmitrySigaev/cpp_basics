#include <iostream>

class A
{
public:
	virtual const char* getName() { return "A"; }
};

class B final : public A // �������� �������� �� ����������� final �����
{
public:
	virtual const char* getName() override { return "B"; }
};

class C : public B // ������ ����������: ������ ����������� ����� final
{
public:
	virtual const char* getName() override { return "C"; }
};

int main()
{
	return 0;
}