
"""

Reference: 

- [Madrigal Database Documentation](https://cedar.openmadrigal.org/docs/name/madContents.html)
- [madrigalWeb](https://github.com/MITHaystack/madrigalWeb): a python API to access the Madrigal database
"""
module MadrigalWeb

# Write your package code here.
using PythonCall
using PythonCall: pynew
using Dates
import Base: basename

const madrigalWeb = pynew()
const MadrigalData = pynew()
const Default_url = "https://cedar.openmadrigal.org"
const Default_server = Ref{Py}()
const User_name = Ref("MadrigalWeb.jl")
const User_email = Ref("")
const User_affiliation = Ref("MadrigalWeb.jl")

export MadrigalData
export download_files
export getExperiments, get_exp_filesExperimentFiles
export filter_by_kindat

function __init__()
    PythonCall.pycopy!(madrigalWeb, pyimport("madrigalWeb.madrigalWeb"))
    PythonCall.pycopy!(MadrigalData, pyimport("madrigalWeb.madrigalWeb").MadrigalData)
end

get_url(server=Default_server[]) = pyconvert(String, server.cgiurl)

function set_default_server(url="https://cedar.openmadrigal.org/")
    get_url(Default_server[]) == url && return Default_server[]
    Default_server[] = MadrigalData(url)
end

function set_default_user(name, email, affiliation=nothing)
    User_name[] = name
    User_email[] = email
    isnothing(affiliation) || (User_affiliation[] = affiliation)
end

abstract type PyObj end


"""
A class that encapsulates information about a Madrigal Experiment.

Similar to the `MadrigalExperiment` class in the madrigalWeb python module.
"""
struct Experiment <: PyObj
    py::Py
end

struct ExperimentFile <: PyObj
    py::Py
end

function Base.getproperty(var::T, s::Symbol) where T<:PyObj
    s in fieldnames(T) ? getfield(var, s) : getproperty(var.py, s)
end

_compat(x::DateTime) = x
_compat(x::String) = DateTime(x)

get_kindat(exp) = pyconvert(Int, exp.kindat)
get_kindatdesc(exp) = pyconvert(String, exp.kindatdesc)
Base.basename(exp::ExperimentFile) = basename(pyconvert(String, exp.name))

function download_files(inst, kindat, t0, t1; dir="./data", server=Default_server[], user_fullname=User_name[], user_email=User_email[], user_affiliation=User_affiliation[])
    t0, t1 = _compat(t0), _compat(t1)
    exps = getExperiments(server, inst, t0, t1)
    files = mapreduce(vcat, exps) do exp
        get_exp_files(server, exp)
    end
    files = filter_by_kindat(files, kindat)
    mkpath(dir)
    map(files) do file
        path = joinpath(dir, basename(file))
        isfile(path) ? path : downloadFile(server, file, path, user_fullname, user_email, user_affiliation)
    end
end

# MadrigalData Class
function downloadFile(server, filename, destination, name=User_name[], email=User_email[], affiliation=User_affiliation[], format="hdf5")
    server.downloadFile(filename, destination, name, email, affiliation, format)
end

function downloadFile(server, expFile::ExperimentFile, destination, name=User_name[], email=User_email[], affiliation=User_affiliation[], format="hdf5")
    server.downloadFile(expFile.name, destination, name, email, affiliation, format)
end

function getExperiments(server, code, startyear, startmonth, startday, starthour, startmin, startsec, endyear, endmonth, endday, endhour, endmin, endsec)
    server.getExperiments(code, startyear, startmonth, startday, starthour, startmin, startsec, endyear, endmonth, endday, endhour, endmin, endsec)
end

function getExperiments(server, code, t0, t1)
    t0 = Dates.DateTime(t0)
    t1 = Dates.DateTime(t1)
    startyear, startmonth, startday, starthour, startmin, startsec = Dates.year(t0), Dates.month(t0), Dates.day(t0), Dates.hour(t0), Dates.minute(t0), Dates.second(t0)
    endyear, endmonth, endday, endhour, endmin, endsec = Dates.year(t1), Dates.month(t1), Dates.day(t1), Dates.hour(t1), Dates.minute(t1), Dates.second(t1)
    server.getExperiments(code, startyear, startmonth, startday, starthour, startmin, startsec, endyear, endmonth, endday, endhour, endmin, endsec)
end


function get_exp_files(server, id::Integer; getNonDefault=false)
    map(ExperimentFile, server.getExperimentFiles(id, getNonDefault))
end

get_exp_files(server, exp; kw...) = map(ExperimentFile, server.getExperimentFiles(exp.id; kw...))

"""
    filter_by_kindat(expFileList, kindat)

Returns a subset of the experiment files in expFileList whose kindat is found in kindat argument.
"""
function filter_by_kindat(expFileList, kindat)
    return filter(expFileList) do expFile
        get_kindat(expFile) in kindat
    end
end
end
