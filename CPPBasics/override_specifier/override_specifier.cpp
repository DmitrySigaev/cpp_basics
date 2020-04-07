#include <iostream>

class A
{
public:
	virtual const char* getName1(int x) { return "A"; }
	virtual const char* getName2(int x) { return "A"; }
};

class B : public A
{
public:
	virtual const char* getName1(short int x) { return "B"; } // ��� ��������� short int
	virtual const char* getName2(int x) const { return "B"; } // ����� �������� const
};

int main()
{
	B b;
	A& rParent = b;
	std::cout << rParent.getName1(1) << '\n';
	std::cout << rParent.getName2(2) << '\n';

	return 0;
}