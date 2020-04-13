#include <iostream>
#pragma pack(1)

struct Empty {};
struct EmptyVirt { virtual ~EmptyVirt() {} };
struct NotEmpty { int m_i; };
struct NotEmptyVirt
{
    virtual ~NotEmptyVirt() {}
    int m_i;
};
struct NotEmptyNonVirt
{
    void foo() const {}
    int m_i;
};


class Parent
{
public:
	// Этот метод getThis() возвращает указатель на класс Parent
	virtual Parent* getThis() { std::cout << "called Parent::getThis()\n"; return this; }
	virtual void printType() { std::cout << "returned a Parent\n"; }
};

class Child : public Parent
{
public:
	// Обычно, типы возврата переопределений и виртуальных функций родительского класса должны совпадать
	// Однако, поскольку Child наследует класс Parent, то следующий метод может возвращать Child* вместо Parent*
	virtual Child* getThis() { std::cout << "called Child::getThis()\n";  return this; }
	void printType() { std::cout << "returned a Child\n"; }
};

int main()
{
	Child ch;
	Parent* p = &ch;
	ch.getThis()->printType(); // вызывается Child::getThis(), возвращается Child*, вызывается Child::printType
	p->getThis()->printType(); // вызывается Child::getThis(), возвращается Parent*, вызывается Parent::printType
	auto apthis = p->getThis();
	apthis->printType();
	((Child *)apthis)->printType();
	p->printType();

    std::cout << sizeof(Empty) << std::endl;
    std::cout << sizeof(EmptyVirt) << std::endl;
    std::cout << sizeof(NotEmpty) << std::endl;
    std::cout << sizeof(NotEmptyVirt) << std::endl;
    std::cout << sizeof(NotEmptyNonVirt) << std::endl;
}