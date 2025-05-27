# MadrigalWeb.jl

[![Build Status](https://github.com/Beforerr/madrigalWeb.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Beforerr/madrigalWeb.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia API to access the Madrigal database.


```julia
using MadrigalWeb
using Dates

MadrigalWeb.set_default_server("https://cedar.openmadrigal.org")
MadrigalWeb.set_default_user("xxx", "xxx@xxx.com", "xxx")

download_files(8100, 10216, "2000-01-01", "2000-01-02")
```