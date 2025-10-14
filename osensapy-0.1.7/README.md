# OsensaPy

Python module for interfacing with the OSENSA transmitters.

## History

_Release 0.1.7 [2019-08-08]:_

* Adjust constructor logic to apply baudrate, parity and timeout settings on port initialization

_Release 0.1.6 [2019-08-08]:_

* Add parameter to allow exclusive port access

_Release 0.1.5 [2019-06-19]:_

* Add method to close transmitter serial port
* Add method to check if serial port is open

_Release 0.1.4 [2018-08-15]:_

* Add logic to flush buffers on read/write

_Release 0.1.3 [2018-04-19]:_

* Add init parameter to control closing serial port between calls

_Release 0.1.2 [2018-03-29]:_

* Fix bug with detecting serial ports in Linux

_Release 0.1.1 [2018-03-02]:_

* Fix issue with internal library import

_Release 0.1.0 [2018-03-01]:_

* Initial release
