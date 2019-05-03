using Geodesics
using SAC
using PlotlyJS

NDKFILE="/Users/zma/remove_noise.dir/cmt2012.ndk"
SACDIR="/Users/zma/remove_noise.dir/test_cascadia.dir/remove_noise_after.dir/"
LOCATIONFILE="/Users/zma/remove_noise.dir/test_cascadia.dir/locations.txt"
REFVEL=4.0

""" parse the ndk file """
function parsendk(filename)
    f=open(filename)
    evtimes=Dict{String,Tuple{Float32,Float32}}()
    while !eof(f)
        tmp=split(readline(f)," ",keepempty=false)
        tmp[2]=replace(tmp[2],"/"=>"-")
        evtimes[string(tmp[2],'T',tmp[3])]=(parse(Float32,tmp[4]),parse(Float32,tmp[5]))
        foreach(x->readline(f),1:4)
    end
    close(f)
    evtimes
end
evtimes=parsendk(NDKFILE)


""" parse the locations file """
function parselocation(filename)
    locs=Dict{String,Tuple{Float32,Float32}}()
    f=open(filename)
    while !eof(f)
        tmp=split(readline(f))
        locs[tmp[2]]=(parse(Float32,tmp[4]),parse(Float32,tmp[5]))
    end
    locs
end
locs=parselocation(LOCATIONFILE)

""" 
compiles the groups for the next commands
each group has the filenames,predicted arrival times
filenames are also sorted by their arrival times 
"""
function makegroups(evtimes,locs,sacdir)
    filenames=filter(x->x[end-2:end]=="SAC",readdir(sacdir))
    groups=Array{Array{Tuple{String,Float32},1},1}()
    
    for key in keys(evtimes)
        println(key)
        lat0,lon0=evtimes[key]
        dists=Array{Float32,1}()
        files=filter(x->occursin(key,x),filenames)
        if length(files)>0
            for filename in files
                nsta=split(filename,'.')[2]
                lat1,lon1=locs[nsta]
                push!(dists,Geodesics.surface_distance(lon0,lat0,lon1,lat1,6371.0)/111.1949)
            end
            idx=sortperm(dists)
            push!(groups,collect(zip(files[idx],dists[idx])))
        end
    end
    groups
end
groups=makegroups(evtimes,locs,SACDIR);

function maketraces(group)
    x=1:7200
    traces=Array{GenericTrace{Dict{Symbol,Any}},1}()
    fac=(group[end][2]-group[1][2])/length(group)
    status.ymin=999
    status.ymax=-999
    for stuff in group
        file,dist=stuff
        tr=SAC.read(SACDIR*file)
        SAC.bp!(tr,0.01,0.1)
        y=(tr.t)[1:7200]
        y=y/maximum(y).*fac.+dist
        if (minimum(y)<status.ymin) status.ymin=minimum(y) end
        if (maximum(y)>status.ymax) status.ymax=maximum(y) end
        push!(traces,scatter(;x=x,y=y,name=file,marker_color="black"))
        push!(traces,scatter(;x=[dist*111.1949/REFVEL],y=[0.5*fac+dist],marker_symbol="triangle-down",
                marker_color="red",mode="markers"))
    end
    traces
end

using PlotlyJS
using Interact
using Blink

const strlen=20   #for the length of buttons

# any better way to save these things????
mutable struct Status
    clicked_points::Array{Tuple{Float32,Float32},1}
    ymin::Float32
    ymax::Float32
end
const status=Status(Array{Tuple{Float32,Float32},1}(),-10.0,10.0)



function genbuttons(contents,callbacks)
    buttons=Array{Widget{:button,Int64},1}()
    for (content,callback) in zip(contents,callbacks)
        if (length(content)>strlen) content=content[1:strlen] end
        nl=div(strlen-length(content),2)
        nr=strlen-nl-length(content)
        push!(buttons,button('_'^nl*content*'_'^nr))
        on(callback,buttons[end])
    end
    buttons
end

# prepare the whole screen
function init_screen()
    layout=Layout(;hovermode="closest",showlegend=false,width=900,height=600)
    win=Window()
    return win,layout
end

global (win,layout)=init_screen();
global p1=plot(scatter(;x=0,y=0),layout);


# here's the part that actually makes the buttons
contents=Array{String,1}()
callbacks=Array{Any,1}()

# button for next
push!(contents,"next")
function next(nclick)
    traces=maketraces(groups[nclick])
    react!(p1,traces,layout)
end
push!(callbacks,next)

# button for select window
push!(contents,"select window")
function selectwindow(ndum)
    relayout!(p1,Dict(:title=>"click on the trace to select window"))
    on(getcursor,p1["click"])
    on(addvertline,p1["click"])
end
function getcursor(data::Dict{String,Any})
    push!(status.clicked_points,
    (Float32(data["points"][1]["x"]),Float32(data["points"][1]["y"])))
end
function addvertline(data)
    if length(status.clicked_points)>2
        deletetraces!(p1,length(p1.plot.data)-1)
        popfirst!(status.clicked_points)
    end
    addtraces!(p1,scatter(;x=[status.clicked_points[end][1],status.clicked_points[end][1]],
        y=[status.ymin,status.ymax],marker_color="red"))
end
push!(callbacks,selectwindow)


# button for confirm and output xpos
push!(contents,"confirm window")
function confirmwindow(ndum)
    if (length(status.clicked_points) != 2)
        relayout!(p1,Dict(:title=>"something not right"))
    else
        relayout!(p1,Dict(:title=>"writing out $(status.clicked_points)"))
        deletetraces!(p1,length(p1.plot.data))
        deletetraces!(p1,length(p1.plot.data))
        pop!(status.clicked_points)
        pop!(status.clicked_points)
        off(p1["click"],getcursor)
        off(p1["click"],addvertline)
    end
end
push!(callbacks,confirmwindow)

# button for delete trace
push!(contents,"delete trace")
function deletetrace(ndum)
    relayout!(p1,Dict(:title=>"click on the trace to delete"))
    on(gettrace,p1["click"])
end
function gettrace(data)
    deletetraces!(p1,data["points"][1]["curveNumber"]+1)   #remember the plus 1
    off(p1["click"],gettrace)
end
push!(callbacks,deletetrace)

buttons=genbuttons(contents,callbacks)

body!(win,vbox(hbox(buttons),p1))
