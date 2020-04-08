#include <iostream>

class A
{
public:
	virtual const char* getName() { return "A"; }
};

class B : public A
{
public:
	// �������� final � �����? ��� ��������, ��� ����� �������������� ��� ������
	virtual const char* getName() override final { return "B"; } // �� ������, ��������������� A::getName()
};

class C : public B
{
public:
	virtual const char* getName() override { return "C"; } // ������ ����������: ��������������� ������ B::getName(), ������� �������� final
};

int main()
{
	return 0;
}