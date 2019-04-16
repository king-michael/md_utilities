# coding=utf-8
import socket


def run_server(fun=print, host='', port=45000, buffer=4096):
    """
    Runs a server and listen to a port.

    Parameters
    ----------
    fun : Callable, optional
        Function to invoke after every recieved data.
        Default is ``print``.
    host : str, optional
        Host address. Use an empty string for localhost. Default is ``''``.
    port : int, optional
        Port to listen to. Default is ``45000``.
    buffer : int, optional
        Buffer size for the stream. Default is ``4096``.

    Returns
    -------
    None
        Runs for ever or till ``KeyboardInterrupt``
    """
    server = socket.socket()
    server.bind((host, port))
    server.listen(1)

    print("Listening on port %d" % port)

    while True:
        try:
            sock, addr = server.accept()
            print("Connection from", sock.getpeername())
            while 1:
                data = sock.recv(buffer)
                # Check if still alive
                if len(data) == 0:
                    break
                # Ignore new lines
                req = data.strip()
                if len(req) == 0:
                    continue
                # Print the request
                #print('Received <--- %s' % req)

                # Do something with it
                #resp = "Hello TCL, this is your response: %s\n" % req.decode('hex')
                #print('Sent     ---> %s' % resp)
                #sock.sendall(resp.encode())

                fun(req)

        except socket.error as ex:
            print('%s' % ex)
            pass
        except KeyboardInterrupt:
            sock.close()
            break


def split_byte_stream(func, conversion=bytes):
    """
    Wrapper to handle a byte stream.
    Splits the byte stream to recieve the different words.
    Invokes ``func`` for every word in the stream.

    Parameters
    ----------
    func : Callable
        Some function
    conversion : Callable
        Conversion function. Default is ``bytes``.

    Returns
    -------
    wrapped_function : Callable
    """
    def wrapper(func):
        def wrapped_func(byte_stream):
            for i in byte_stream.split():
                func(conversion(i))
        return wrapped_func
    return wrapper


if __name__ == '__main__':
    def dummy_function(astr):
        print("vmd send:", astr)

    run_server(dummy_function)
