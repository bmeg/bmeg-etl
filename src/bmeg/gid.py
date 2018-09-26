class GID(str):
    def __new__(cls, content):
        return str.__new__(cls, content)
