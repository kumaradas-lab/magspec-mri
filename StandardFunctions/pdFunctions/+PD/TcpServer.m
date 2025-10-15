classdef TcpServer < handle
  %% TCP (socket) server
  %
  %   ts = basic.TcpServer(options)
  %
  % Start a TCP server with methods for sending and receiving (EOL terminated)
  % strings.
  % The construction of an object of this class blocks the execution until a
  % client has connected or the timeout is exceeded (see below).
  %
  % INPUT:
  %
  %   options
  %         Optional. Structure with settings for the TCP server. If fields are
  %         omitted or empty, default values are used:
  %
  %     port
  %           Port on the localhost for the TCP server. If this port cannot be
  %           opened an error is thrown. (Default: 4444)
  %
  %     timeout
  %           Timeout for the client connection in ms. If no client connects
  %           before the timeout is reached, an error is thrown. If this is set
  %           to 0, the server waits for a client indefinitely. (Default: 10000)
  %
  % ----------------------------------------------------------------------------
  % (C) Copyright 2020-2023 Pure Devices GmbH, Wuerzburg, Germany
  % www.pure-devices.com
  % ----------------------------------------------------------------------------


  properties (SetAccess = immutable)
    port
    timeout
  end


  properties (Access = protected)
    serverSocket
    clientSocket
    inputStream
    outputStream
  end


  methods

    function this = TcpServer(options)
      %% Open TCP server and wait for a client to connect

      if nargin < 1 || isempty(options), options = struct(); end

      if isemptyfield(options, 'port'), options.port = 4444; end  % port for the TCP server
      if isemptyfield(options, 'timeout'), options.timeout = 10000; end  % timeout for client connection in ms

      if ~isnumeric(options.port) || mod(options.port, 1) ~= 0 || ...
          options.port < 1025 || options.port > 65535
        error('PD:TcpServer:InvalidPort', 'Invalid port');
      end
      if ~isnumeric(options.timeout) || mod(options.timeout, 1) ~= 0 || ...
          options.timeout < 0
        error('PD:TcpServer:InvalidTimeout', 'Invalid timeout');
      end

      this.port = options.port;
      this.timeout = options.timeout;

      this.serverSocket = java.net.ServerSocket(this.port);

      this.serverSocket.setSoTimeout(this.timeout);
      % This command blocks until a client has connected or timeout is exceeded
      this.clientSocket = this.serverSocket.accept();
      this.inputStream = this.clientSocket.getInputStream();
      this.outputStream = this.clientSocket.getOutputStream();
    end


    function delete(this)
      %% Close server on destruction of object

      if ~isempty(this.serverSocket) && ~this.serverSocket.isClosed()
        this.serverSocket.close();
      end
    end


    function str = readlNonBlocking(this)
      %% Read available bytes sent by client as string if available
      %
      % All available bytes or the bytes until an EOL character (whichever comes
      % first) are returned as string.

      str = '';
      while (this.inputStream.available)
        % cast received byte to char
        str(end+1) = char(this.inputStream.read());
        if str(end) == sprintf('\n')
          break;
        end
      end
    end


    function str = receivelNonBlocking(this)
      %% Read string from client until next EOL if available
      %
      % Only if a message is available, this function blocks execution until an
      % EOL character is received. It is non-blocking if no message is
      % available at all.

      str = '';
      while (true)
        str = [str, this.readlNonBlocking()];
        if isempty(str) || str(end) == sprintf('\n')
          break;
        end
        java.lang.Thread.sleep(10);
      end
    end


    function str = readl(this)
      %% Read bytes sent by client as string if available
      %
      % All available bytes or the bytes until an EOL character (whichever comes
      % first) are returned as string.
      % This function blocks execution until a message is received

      str = char(this.inputStream.read());
      while (this.inputStream.available)
        % cast received byte to char
        str(end+1) = char(this.inputStream.read());
        if str(end) == sprintf('\n')
          break;
        end
      end
    end


    function str = receivel(this)
      %% Read string from client until next EOL if available
      %
      % This function returns when an EOL character is received.
      % It throws an error if the client closed the connection.

      this.clientSocket.setSoTimeout(2);

      str = '';
      while (true)
        byte = 0;
        try %#ok<TRYNC>
          byte = this.inputStream.read();
          str = [str, char(byte)];
        end
        if byte == -1
          this.clientSocket.setSoTimeout(0);
          error('PD:TcpServer:ConnectionClosed', ...
            'The client seems to have closed the connection');
        end
        str = [str, this.readlNonBlocking()];
        if ~isempty(str) && str(end) == sprintf('\n')
          break;
        end
        java.lang.Thread.sleep(10);
      end

      this.clientSocket.setSoTimeout(0);

    end


    function write(this, str)
      %% Send (EOL terminated) string to client

      if isempty(str)
        return;
      end

      if str(end) ~= sprintf('\n')
        % append EOL
        str(end+1) = sprintf('\n');
      end
      fprintf('%s: sending message "%s"...\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'), lower(str(1:end-1))),
      % send as bytes
      this.outputStream.write(uint8(str));
    end

  end

end
