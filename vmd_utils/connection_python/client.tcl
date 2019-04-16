# VMD Script
# Details
# -------
# sends frame from vmd to python

namespace eval ::connection_python:: {
  variable socket_host "127.0.0.1"
  variable socket_port 45000
  variable sock ""
  variable trace 0

  proc help {} {
    variable socket_host
    variable socket_port
    variable sock
    puts "Proc to connect to python and send frames to it"
    puts "\nConfiguration:\n"
    puts "  host : ${socket_host}"
    puts "  port : ${socket_port}"
    puts "\nCommands:\n"
    puts "::connection_python::config ?host? ?port?"
    puts "  Change host and port"
    puts "::connection_python::connect"
    puts "   Connects to python socket"
    puts "::connection_python::disconnect"
    puts "   Disconnects from the python socket"
    puts "::connection_python::send_frame frame"
    puts "   Sends frame to python"
    puts "::connection_python::start ?molid?"
    puts "   Trace vmd_frame and send it to python"
    puts "::connection_python::stop ?molid?"
    puts "   Stop tracing vmd_frame"
    puts "::connection_python::restart ?molid?"
    puts "   Restart tracing vmd_frame"
  }

  proc config {{host "127.0.0.1"} {port 45000}} {
    # set some globals
    variable socket_host ${host}
    variable socket_port ${port}
  }
  proc send_frame {frame} {
    variable sock
    puts $sock "${frame}"
  }

  proc send_frame_trace {name index op} {
    # name == vmd_frame
    # index == molecule id of the newly changed frame
    # op == w
    set frame [molinfo $index get frame]
    ::connection_python::send_frame ${frame}
  }

  proc connect {} {
    variable socket_host
    variable socket_port
    variable sock
    if {[llength ${sock}]} {
      puts "connection_python) already connected to : ${socket_host} ${socket_port}"
      return
    }
    puts "connection_python) connect to : ${socket_host} ${socket_port}"
    # Connect to server
    set sock [socket $socket_host $socket_port]
    # Disable line buffering
    fconfigure $sock -buffering none
  }

  proc disconnect {} {
    # close socket
    variable sock
    if {![string compare ${sock} ""]} {
      puts "connection_python) already diconnected"
      return}
    close $sock
    set sock ""
    puts "connection_python) disconnect from socket"
  }

  proc start {{molid top}} {
    variable trace
    if {${trace}} {::connection_python::stop}
    ::connection_python::connect

    # connect frame change
    global vmd_frame
    if {! [string compare $molid top]} {set molid [molinfo top]}
    trace variable vmd_frame($molid) w ::connection_python::send_frame_trace
    set trace 1; puts "connection_python) start trace vmd_frame"
  }

  proc stop {{molid top}} {
    # disconnect frame change
    global vmd_frame
    if {! [string compare $molid top]} {set molid [molinfo top]}
    trace vdelete vmd_frame($molid) w ::connection_python::send_frame_trace
    variable trace 0; puts "connection_python) stop trace vmd_frame"
    ::connection_python::disconnect
  }

  proc restart {} {
    ::connection_python::stop
    ::connection_python::start
  }

}
proc connection_python {{cmd help} args} {
  if {[llength $args]} {
    ::connection_python::${cmd} {*}$args
  } else {
    ::connection_python::${cmd}
  }
}
