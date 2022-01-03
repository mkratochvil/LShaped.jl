filestring = string(@__FILE__)
n = length(filestring)
if filestring[n-4] == 'e'
    arrayid = parse(Int64,filestring[n-3])
    println(arrayid)
    println(typeof(arrayid))
elseif filestring[n-5] == 'e'
    arrayid = parse(Int64,filestring[n-4:n-3])
    println(arrayid)
    println(typeof(arrayid))
elseif filestring[n-6] == 'e'
    arrayid = parse(Int64,filestring[n-5:n-3])
    println(arrayid)
    println(typeof(arrayid))
elseif filestring[n-7] == 'e'
    arrayid = parse(Int64,filestring[n-6:n-3])
    println(arrayid)
    println(typeof(arrayid))
end