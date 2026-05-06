using DMSP
using Dates
id, trange = 16, [DateTime("2022-07-02T05:00"), DateTime("2022-07-02T05:10")]

DMSP.flux(trange, id)