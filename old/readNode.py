class ReadNode():
    '''
    Represents a node in a linked list of reads
    '''

    def __init__(self, pos, pair=None, next=None):
        self.pos = pos
        self.pair = pair
        self.next = next

    def addRead(self, node):
        if self.pos > node.pos:
            # Add new node at beginning
            node.next = self
            return node
        else:
            curr = self
            while curr.next and curr.next.pos <= node.pos:
                curr = curr.next

            node.next = curr.next
            curr.next = node

            return self

    def index(self, id=0):
        #print('  Indexing #%d' % id)
        self.id = id

        curr = self
        while curr:
            curr.id = id
            id += 1
            curr = curr.next

    def toString(self):
        if self.next:
            return str(self.pos) + ' --> ' + self.next.toString()
        else:
            return str(self.pos)
