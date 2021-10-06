/**
 * Source Code first verified at https://etherscan.io on Thursday, March 21, 2019
 (UTC) */

pragma solidity 0.4.24;

library SafeMath {

    function mul(uint256 a, uint256 b) internal pure returns (uint256 c) {
        if (a == 0) {
            return 0;
        }
        c = a * b;
        assert(c / a == b);
        return c;
    }

    function div(uint256 a, uint256 b) internal pure returns (uint256) {
        return a / b;
    }

    function sub(uint256 a, uint256 b) internal pure returns (uint256) {
        assert(b <= a);
        return a - b;
    }

    function add(uint256 a, uint256 b) internal pure returns (uint256 c) {
        c = a + b;
        assert(c >= a);
        return c;
    }
}

contract Ownable {
    address public owner;
    event OwnershipTransferred(address indexed previousOwner, address indexed newOwner);

    constructor() public {
        owner = msg.sender;
    }

    modifier onlyOwner() {
        require(msg.sender == owner);
        _;
    }

    function transferOwnership(address newOwner) public onlyOwner {
        require(newOwner != address(0));
        emit OwnershipTransferred(owner, newOwner);
        owner = newOwner;
    }

}

contract ERC20Basic {
    function totalSupply() public view returns (uint256);
    function balanceOf(address who) public view returns (uint256);
    function transfer(address to, uint256 value) public returns (bool);
    event Transfer(address indexed from, address indexed to, uint256 value);
}

contract ERC20 is ERC20Basic {
    function allowance(address owner, address spender) public view returns (uint256);
    function transferFrom(address from, address to, uint256 value) public returns (bool);
    function approve(address spender, uint256 value) public returns (bool);
    event Approval(address indexed owner, address indexed spender, uint256 value);
}

contract BasicToken is ERC20Basic {
    using SafeMath for uint256;
    
    mapping(address => uint256) balances;
    
    uint256 totalSupply_;
    
    function totalSupply() public view returns (uint256) {
        return totalSupply_;
    }

    function transfer(address _to, uint256 _value) public returns (bool) {
        require(_to != address(0));
        require(_value <= balances[msg.sender]);

        balances[msg.sender] = balances[msg.sender].sub(_value);
        balances[_to] = balances[_to].add(_value);
        emit Transfer(msg.sender, _to, _value);
        return true;
    }

    function balanceOf(address _owner) public view returns (uint256) {
        return balances[_owner];
    }

}

contract StandardToken is ERC20, BasicToken {

    mapping (address => mapping (address => uint256)) internal allowed;

    function transferFrom(address _from, address _to, uint256 _value) public returns (bool) {
        require(_to != address(0));
        require(_value <= balances[_from]);
        require(_value <= allowed[_from][msg.sender]);

        balances[_from] = balances[_from].sub(_value);
        balances[_to] = balances[_to].add(_value);
        allowed[_from][msg.sender] = allowed[_from][msg.sender].sub(_value);
        emit Transfer(_from, _to, _value);
        return true;
    }

    function approve(address _spender, uint256 _value) public returns (bool) {
        allowed[msg.sender][_spender] = _value;
        emit Approval(msg.sender, _spender, _value);
        return true;
    }

    function allowance(address _owner, address _spender) public view returns (uint256) {
        return allowed[_owner][_spender];
    }

    function increaseApproval(address _spender, uint _addedValue) public returns (bool) {
        allowed[msg.sender][_spender] = allowed[msg.sender][_spender].add(_addedValue);
        emit Approval(msg.sender, _spender, allowed[msg.sender][_spender]);
        return true;
    }

    function decreaseApproval(address _spender, uint _subtractedValue) public returns (bool) {
        uint oldValue = allowed[msg.sender][_spender];
        if (_subtractedValue > oldValue) {
            allowed[msg.sender][_spender] = 0;
        } else {
            allowed[msg.sender][_spender] = oldValue.sub(_subtractedValue);
        }
        emit Approval(msg.sender, _spender, allowed[msg.sender][_spender]);
        return true;
    }

    function _burn(address account, uint256 amount) internal {
        require(account != 0);
        require(amount <= balances[account]);

        totalSupply_ = totalSupply().sub(amount);
        balances[account] = balances[account].sub(amount);
        emit Transfer(account, address(0), amount);
    }

    function _burnFrom(address account, uint256 amount) internal {
        require(amount <= allowed[account][msg.sender]);
        allowed[account][msg.sender] = allowed[account][msg.sender].sub(amount);
        _burn(account, amount);
    }

}


contract BurnableToken is StandardToken, Ownable {

  function burn(uint256 value) public {
    _burn(msg.sender, value);
  }

  function burnFrom(address from, uint256 value) public {
    _burnFrom(from, value);
  }

  function _burn(address who, uint256 value) internal {
    super._burn(who, value);
  }
}

contract initialSupplyToken is StandardToken, Ownable {
    using SafeMath for uint256;

    function initialSupply(uint256 _amount) public onlyOwner returns (bool) {
        totalSupply_ = totalSupply_.add(_amount);
        balances[msg.sender] = balances[msg.sender].add(_amount);
        return true;
    }

}
contract BidaCoin is  initialSupplyToken,BurnableToken  {
    using SafeMath for uint256;

    string public name = "BidaCoin";
    string public symbol = "BIDA";
    uint8 constant public decimals = 18;
    
    constructor () public {
        initialSupply(300000000 * (10 ** uint256(decimals)));
    }

    struct LockParams {
        uint256 TIME;
        uint256 AMOUNT;
    }

    mapping(address => LockParams[]) private holdAmounts;
    address[] private holdAmountAccounts;

    function isValidAddress(address _address) public view returns (bool) {
        return (_address != 0x0 && _address != address(0) && _address != 0 && _address != address(this));
    }

    modifier validAddress(address _address) {
        require(isValidAddress(_address));
        _;
    }


    function transfer(address _to, uint256 _value) public validAddress(_to) returns (bool) {
        require(checkAvailableAmount(msg.sender, _value));

        return super.transfer(_to, _value);
    }

    function transferFrom(address _from, address _to, uint256 _value) public validAddress(_to) returns (bool) {
        require(checkAvailableAmount(_from, _value));

        return super.transferFrom(_from, _to, _value);
    }

    function approve(address _spender, uint256 _value) public returns (bool) {
        return super.approve(_spender, _value);
    }

    function setHoldAmount(address _address, uint256 _amount, uint256 _time) public onlyOwner {
        _setHold(_address, _amount, _time);
    }

    function _setHold(address _address, uint256 _amount, uint256 _time) internal {
        LockParams memory lockdata;
        lockdata.TIME = _time;
        lockdata.AMOUNT = _amount;
        holdAmounts[_address].push(lockdata);
        holdAmountAccounts.push(_address) - 1;
    }

    function getTotalHoldAmount(address _address) public view returns(uint256) {
        uint256 totalHold = 0;
        LockParams[] storage locks = holdAmounts[_address];
        for (uint i = 0; i < locks.length; i++) {
            if (locks[i].TIME >= now) {
                totalHold = totalHold.add(locks[i].AMOUNT);
            }
        }
        return totalHold;
    }
    
    function getNow() public view returns(uint256) {
        return now;
    }

    function getAvailableBalance(address _address) public view returns(uint256) {
        return balanceOf(_address).sub(getTotalHoldAmount(_address));
    }

    function checkAvailableAmount(address _address, uint256 _amount) public view returns (bool) {
        return _amount <= getAvailableBalance(_address);
    }

    function removeHoldByAddress(address _address) public onlyOwner {
        delete holdAmounts[_address];
    }

    function removeHoldByAddressIndex(address _address, uint256 _index) public onlyOwner {
        LockParams[] memory nowLocks = holdAmounts[_address];
        removeHoldByAddress(_address);
        for (uint i = 0; i < nowLocks.length; i++) {
            if (i != _index) {
                 setHoldAmount(_address,nowLocks[i].AMOUNT,nowLocks[i].TIME);
            }
        }
    }


    function  changeHoldByAddressIndex(address _address, uint256 _index, uint256 _amount, uint256 _time) public onlyOwner {
        holdAmounts[_address][_index].TIME = _time;
        holdAmounts[_address][_index].AMOUNT = _amount;
    }


    function getHoldAmountAccounts() public view returns (address[]) {
        return holdAmountAccounts;
    }

    function countHoldAmount(address _address) public view returns (uint256) {
        require(_address != 0x0 && _address != address(0));
        return holdAmounts[_address].length;
    }

    function getHoldAmount(address _address, uint256 _idx) public view  returns (uint256, uint256) {
        require(_address != 0x0);
        require(holdAmounts[_address].length>0);
        return (holdAmounts[_address][_idx].TIME, holdAmounts[_address][_idx].AMOUNT);
    }

    function burnFrom(address from, uint256 value) public {
        require(checkAvailableAmount(from, value));
        super.burnFrom(from, value);
    }

    function burn(uint256 value) public {
        require(checkAvailableAmount(msg.sender, value));
        super.burn(value);
    }

}