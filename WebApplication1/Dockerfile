FROM microsoft/aspnetcore:2.0 AS base
WORKDIR /app
EXPOSE 80

FROM microsoft/aspnetcore-build:2.0 AS build
WORKDIR /src
COPY CiCd.sln ./
COPY WebApplication1/WebApplication1.csproj WebApplication1/
COPY WebApplication1.Tests/WebApplication1.Tests.csproj WebApplication1.Tests/
RUN dotnet restore -nowarn:msb3202,nu1503
COPY . .
WORKDIR /src/WebApplication1
RUN dotnet build -c Release -o /app

FROM build as test
WORKDIR /src/WebApplication1.Tests

FROM build AS publish
WORKDIR /src/WebApplication1
RUN dotnet publish -c Release -o /app

FROM base AS final
WORKDIR /app
COPY --from=publish /app .
ENTRYPOINT ["dotnet", "WebApplication1.dll"]
