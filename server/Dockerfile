# На всякий случай напоминаю, что Dockerfile всегда должен начинаться с импорта базового образа.
# Для этого используется ключевое слово 'FROM'.
# Здесь нам нужно импортировать образ python (с DockerHub).
# В результате мы, в качестве имени образа, указываем 'python', а в качестве версии - 'latest'.

FROM python:latest
#FROM sixcomp/sixcom:latest
ENV DEBIAN_FRONTEND noninteractive


#we rollback when it build successful
#RUN apt-get update && \
#       apt-get install -y --no-install-recommends apt-transport-https ca-certificates gnupg software-properties-common wget

RUN apt-get update && \
    apt-get install -y --no-install-recommends zip unzip git build-essential  ca-certificates \
     python3-dev make curl tar software-properties-common qtbase5-dev libgsl-dev  libgtest-dev libboost-all-dev

#RUN  wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null

#RUN  apt-add-repository 'deb https://apt.kitware.com/ubuntu/ focal main'

#RUN  apt-get update

#RUN  apt-get install -y --no-install-recommends kitware-archive-keyring
#RUN  rm /etc/apt/trusted.gpg.d/kitware.gpg

#RUN  apt-get install -y --no-install-recommends cmake

RUN curl -s "https://cmake.org/files/v3.17/cmake-3.17.0-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local
#libx11-dev libmesa-dev libxi-dev libxext-dev

#RUN curl -s http://ftp.us.debian.org/debian/pool/main/g/gcc-9/g++-9_9.3.0-13_amd64.deb | tar --strip-components=1 -xz -C /usr/local
#
   
RUN  echo "deb http://deb.debian.org/debian experimental main" | tee /etc/apt/sources.list.d/docker.list
                                                                          
RUN apt-get update -y && apt-get -t experimental install -y --no-install-recommends g++-9


#RUN cd /opt && git clone https://github.com/microsoft/vcpkg.git vcpkg
#RUN cd /opt/vcpkg && ls .
#RUN cd /opt/vcpkg && ./bootstrap-vcpkg.sh

####RUN cd /opt/vcpkg && ./vcpkg install python3:x64-linux
#RUN cd /opt/vcpkg && ./vcpkg install boost-python:x64-linux
#RUN cd /opt/vcpkg && ./vcpkg install gsl:x64-linux
#RUN cd /opt/vcpkg && ./vcpkg install gtest:x64-linux
#RUN cd /opt/vcpkg && ./vcpkg install qt5[core]:x64-linux

#tests

# Для того чтобы запустить в контейнере код, написанный на Python, нам нужно импортировать файлы 'server.py' и 'index.html'.
# Для того чтобы это сделать, мы используем ключевое слово 'ADD'.
# Первый параметр, 'server.py', представляет собой имя файла, хранящегося на компьютере.
# Второй параметр, '/server/', это путь, по которому нужно разместить указанный файл в образе.
# Здесь мы помещаем файл в папку образа '/server/'.

ADD server.py /server/
ADD index.html /server/

# Здесь мы воспользуемся командой 'WORKDIR', возможно, новой для вас.
# Она позволяет изменить рабочую директорию образа.
# В качестве такой директории, в которой будут выполняться все команды, мы устанавливаем '/server/'.

WORKDIR /server/