from __future__ import absolute_import
import requests
import configparser  # @todo Should check for this in __init__?
import json
from os import environ, path
from getpass import getpass
import base64


def get_patric_token(username, password=None, token_type='rast', timeout=5):
    """ Get an authentication token for Patric web services.

        @param username: User name
        @param password: Password or None to prompt and enter password
        @param token_type: Type of authentication token
        @param timeout:  Number of seconds to wait for response
        @return: Authentication token
    """

    # Prompt for a password if not specified.
    if password is None:
        password = getpass(prompt='{0} password: '.format(token_type))

    # Get an authentication token from the specified web service.
    if token_type == 'rast':
        url = 'http://rast.nmpdr.org/goauth/token?grant_type=client_credentials&client_id={0}'.format(username)
        headers = dict()
        headers['Content-Type'] = 'application/json'
        headers['Authorization'] = 'Basic ' + base64.encodestring(username + ':' + password)[:-1]
        response = requests.get(url, headers=headers, timeout=timeout)
        if response.status_code != requests.codes.OK:
            response.raise_for_status()
        data = json.loads(response.text)
        token = data['access_token']
        user_id = data['client_id']

    else:
        raise ValueError('Token token_type {0} is not valid'.format(token_type))

    # Save the authentication data in config file.
    config_file = path.join(environ['HOME'], '.patric_config')
    config = configparser.ConfigParser()
    config.read(config_file)
    if not config.has_section('authentication'):
        config.add_section('authentication')
    config.set('authentication', 'token', token)
    config.set('authentication', 'user_id', user_id)
    config.write(open(config_file, 'w'))

    return token


def shock_download(url, token):
    """ Download data from a Shock node.

        @param url: URL to Shock node
        @param token: Authentication token for Patric web services
        @return Data from Shock node
    """

    response = requests.get(url + '?download', headers={'AUTHORIZATION': 'OAuth ' + token})
    if response.status_code != requests.codes.OK:
        response.raise_for_status()
    return response.text


class ServerError(Exception):
    """ Exception raised when server returns an error. """

    def __init__(self, message, name='JSONRPCError', code=-32603, data=None, error=None):
        # The name is always 'JSONRPCError and the code is always -32603 because it
        # is hard-coded in the server.
        self.name = name
        self.code = code

        if message is None:
            self.data = data or error or ''
            self.message = self.data
        else:
            # Separate the lines in the message and store as an array of lines.
            self.data = message.split('\n')

            # If the first line has the _ERROR_ delimiter, extract the message.
            if self.data[0].startswith('_ERROR_'):
                self.message = self.data[0][7:-7]

            # If the second line has the _ERROR_ delimiter, extract the message.
            if self.data[0] == 'JSONRPC error:' and len(self.data) > 1 and self.data[1].startswith('_ERROR_'):
                self.message = self.data[1][7:-7]

            # Otherwise, use the first line as the message.
            else:
                self.message = self.data[0]

    def traceback(self):
        output = ''
        for line in self.data:
            output += line + '\n'
        return output

    def __str__(self):
        return self.message


class ObjectNotFoundError(Exception):
    """ Exception raised when an object cannot be retrieved from workspace. """

    def __init__(self, message, data):
        self.message = message
        self.data = data

    def traceback(self):
        output = ''
        for line in self.data:
            output += line + '\n'
        return output

    def __str__(self):
        return self.message


class ObjectTypeError(Exception):
    """ Exception raised when an object is the wrong type for the operation. """

    def __init__(self, message, data):
        self.message = message
        self.data = data

    def traceback(self):
        output = ''
        for line in self.data:
            output += line + '\n'
        return output

    def __str__(self):
        return 'ObjectTypeError: ' + self.message


class DuplicateGapfillError(Exception):
    """ Exception raised when a gap fill solution is already available. """
    pass


class JobError(Exception):
    """ Exception raised when a job submitted to app service ended in error. """
    pass


class PatricClient(object):
    """ Client for Patric web services """

    def __init__(self, url, name, token=None):
        """ Initialize object.

            @param url: URL of service endpoint
            @param name: Name of service
            @param token: Authentication token for Patric web services
        """

        self.url = url
        self.name = name

        if token is None:
            # Get the authentication token from the Patric config file.
            config = configparser.ConfigParser()
            config.read(path.join(environ['HOME'], '.patric_config'))
            self.token = config.get('authentication', 'token')
            self.user_id = config.get('authentication', 'user_id')
        else:
            self.token = token
            self.user_id = token.split('|')[0].replace('un=', '')

        # Create the headers for the request to the server.
        self.headers = dict()
        self.headers['AUTHORIZATION'] = self.token

        return

    def call(self, method, params, timeout=1800):
        """ Call a server method and wait for the response.

            @param method: Name of method
            @param params: Dictionary of input parameters for method
            @param timeout: Number of seconds to wait for response
            @return Output of method in JSON format
        """

        # Create the body of the request for the specified method.
        request_data = dict()
        request_data['method'] = self.name + '.' + method
        request_data['params'] = params
        request_data['version'] = '1.1'
        request_data['id'] = '1'

        # Send the request to the server and get back a response.
        response = requests.post(self.url, data=json.dumps(request_data), headers=self.headers, timeout=timeout)

        if response.status_code == requests.codes.server_error:
            if 'content-type' in response.headers and response.headers['content-type'] == 'application/json':
                err = json.loads(response.text)
                if 'error' in err:
                    raise ServerError(**err['error'])
                else:
                    raise ServerError(response.text)
            else:
                raise ServerError(response.text)

        if response.status_code != requests.codes.OK:
            response.raise_for_status()
        return json.loads(response.text)['result'][0]  # Get the output from the method in the response
